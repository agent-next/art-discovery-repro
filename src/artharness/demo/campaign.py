"""Wire a whole mini-campaign: world, briefs, gates, triage, agents, report, scorecard.

This is the paper's campaign at toy scale with the *same* control flow: five gated
stages, worker/supervisor/curator sessions per task, follow-ups through triage, an
editor before the report is filed, git-versioned task records and session accounting.
Only the scale (9 contigs, ~15 tasks) and the model (a rule policy or a 0.6B local
model) differ.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from ..accounting import SessionLedger
from ..config import CampaignConfig
from ..knowledge import KnowledgeBase
from ..orchestrator import STAGES, CampaignReport, Orchestrator
from ..records import RecordStore, TaskStatus
from ..runner.base import ScriptedBackend  # noqa: F401  (re-export for tests)
from .agents import OllamaAgent, OllamaChat, RuleAgent, check_ollama
from .narrator import NarratedBackend, Narrator
from .replay import write_replay
from .tools import FRAGMENT_MAX_AA
from .walkthrough import write_walkthrough
from .world import build_world, read_fasta

BRIEFS = {
    "1_input_assembly": [
        "Stage 1 - input assembly. Inventory the input contigs: how many and how long.\n"
        "Tools: list_contigs"],
    "2_database_sweep": [
        "Stage 2 - database sweep. Search all six reading frames of every contig for open "
        "reading frames that carry the RT catalytic motif.\nTools: find_rt_orfs"],
    "3_rt_classification": [
        "Stage 3 - RT classification. Give every RT candidate a class and set fragments "
        f"(under {FRAGMENT_MAX_AA} aa) aside.\nTools: classify_rt"],
    "4_neighborhood_census": [
        "Stage 4 - neighborhood census. Scan the window upstream of every full-length RT "
        "for a tandem repeat array (paper Methods p.32 step 1, with shuffle control). "
        "Propose one follow-up deep dive per array call.\nTools: scan_arrays"],
    "5_deep_dives": [
        "Stage 5 - deep dives. Make sure every array-positive locus has a deep dive; run "
        "one for any that lacks it.\nTools: deep_dive"],
}
_LOCUS = re.compile(r"\bL\d{2}\b")


@dataclass
class DemoOptions:
    out: Path
    agent: str = "rules"  # rules | ollama
    model: str = "qwen3:0.6b"
    ollama_url: str = "http://localhost:11434"
    stage_revision: bool = True
    concurrency: int = 1
    explain: bool = True
    force: bool = False
    seed: int = 20260929


@dataclass
class Scorecard:
    true_arrays: int = 0
    found_true: int = 0
    false_calls: list[str] = field(default_factory=list)
    missed: list[str] = field(default_factory=list)
    decoys_ok: dict[str, bool] = field(default_factory=dict)

    @property
    def perfect(self) -> bool:
        return not self.false_calls and not self.missed and all(self.decoys_ok.values())


@dataclass
class DemoResult:
    out: Path
    report: CampaignReport
    scorecard: Scorecard | None
    ok: bool
    error: str | None
    filed: Path | None
    walkthrough: Path


def _table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    cols = lines[0].split("\t")
    return [dict(zip(cols, ln.split("\t"), strict=True)) for ln in lines[1:] if ln]


def make_gates(shared: Path, genomes: Path, narrator: Narrator) -> dict[str, Callable]:
    """Scripted completion checks. Each recomputes what it needs from the raw input, so
    it never trusts an agent's word for the thing it is checking."""

    def n_contigs() -> dict[str, int]:
        return {k: len(v) for k, v in read_fasta(genomes).items()}

    def g1() -> str | None:
        rows = _table(shared / "contigs.tsv")
        want = n_contigs()
        return None if len(rows) == len(want) and rows else \
            f"contigs.tsv has {len(rows)} rows, input has {len(want)}"

    def g2() -> str | None:
        rows, lens = _table(shared / "rt_candidates.tsv"), n_contigs()
        bad = [r["locus_id"] for r in rows if r["contig"] not in lens
               or not 0 <= int(r["start"]) < int(r["end"]) <= lens[r["contig"]]]
        return "no RT candidates" if not rows else (f"bad coordinates: {bad}" if bad else None)

    def g3() -> str | None:
        cands = {r["locus_id"] for r in _table(shared / "rt_candidates.tsv")}
        got = {r["locus_id"] for r in _table(shared / "rt_classes.tsv")}
        return None if cands and cands == got else \
            f"classes cover {len(got)} of {len(cands)} candidates"

    def g4() -> str | None:
        full = {r["locus_id"] for r in _table(shared / "rt_classes.tsv") if r["fragment"] == "0"}
        rows = _table(shared / "arrays.tsv")
        got = [r["locus_id"] for r in rows]
        ok_status = all(r["status"] in ("array", "no_array", "not_assessed") for r in rows)
        return None if full and set(got) == full and len(got) == len(full) and ok_status \
            else f"arrays.tsv covers {len(set(got))} of {len(full)} full-length RTs"

    def g5() -> str | None:
        arr = [r["locus_id"] for r in _table(shared / "arrays.tsv") if r["status"] == "array"]
        missing = [a for a in arr if not (shared / "deep_dives" / f"{a}.tsv").exists()]
        return f"no deep dive for {missing}" if missing else None

    checks = dict(zip(STAGES, (g1, g2, g3, g4, g5), strict=True))

    def wrap(stage: str, check: Callable[[], str | None]) -> Callable[[str], bool]:
        def gate(_stage: str) -> bool:
            why = check()
            narrator.event("gate", stage=stage, ok=why is None, why=why or "")
            narrator.emit(f"  GATE {stage}: " + ("PASS" if why is None else f"FAIL - {why}"),
                          "gate")
            return why is None
        return gate

    return {s: wrap(s, c) for s, c in checks.items()}


def make_triage(shared: Path, narrator: Narrator) -> Callable:
    opened: set[str] = set()

    def triage(brief: str, parent) -> tuple[bool, str]:
        m = _LOCUS.search(brief)
        if not m:
            ok, reason = False, "the follow-up names no locus, so there is nothing concrete to run"
        else:
            lid = m.group(0)
            status = {r["locus_id"]: r["status"] for r in _table(shared / "arrays.tsv")}
            if lid in opened:
                ok, reason = False, f"a deep dive on {lid} is already queued"
            elif status.get(lid) != "array":
                ok, reason = False, (f"arrays.tsv calls {lid} '{status.get(lid, 'absent')}'; "
                                     "a deep dive needs an array call")
            else:
                opened.add(lid)
                ok, reason = True, ""
        narrator.event("triage", parent=parent.task_id, brief=brief[:140], ok=ok, why=reason)
        narrator.emit(f"  TRIAGE from {parent.task_id}: "
                      + ("released" if ok else f"REJECTED - {reason}") + f"  [{brief[:60]}]",
                      "triage")
        return ok, reason

    return triage


def build_report(shared: Path, task_ids: list[str]) -> str:
    cands = {r["locus_id"]: r for r in _table(shared / "rt_candidates.tsv")}
    classes = {r["locus_id"]: r for r in _table(shared / "rt_classes.tsv")}
    arrays = {r["locus_id"]: r for r in _table(shared / "arrays.tsv")}
    dives = {p.stem: _table(p)[0] for p in sorted((shared / "deep_dives").glob("*.tsv"))}
    pos = [lid for lid, r in arrays.items() if r["status"] == "array"]
    lines = [
        "# Array-associated RT candidates in the toy world",
        "",
        "Demo campaign; the brief is a reconstruction, not Anthropic's text. "
        f"Tasks: {', '.join(task_ids)}.",
        "",
        "## Result",
        f"{len(pos)} array-positive loci among {len(arrays)} assessed full-length RTs "
        f"({len(cands)} RT candidates in total).",
        "",
        "| locus | contig | class | R | deep dive | conserved block nt | median spacing nt |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for lid in pos:
        d = dives.get(lid, {})
        lines.append(f"| {lid} | {cands[lid]['contig']} | {classes[lid]['rt_class']} | "
                     f"{arrays[lid]['R']} | {d.get('status', 'missing')} | "
                     f"{d.get('block_len', '-')} | {d.get('median_spacing_nt', '-')} |")
    lines += ["", "## Set aside"]
    frag = [lid for lid, c in classes.items() if c["fragment"] == "1"]
    na = [lid for lid, a in arrays.items() if a["status"] == "not_assessed"]
    none = [lid for lid, a in arrays.items() if a["status"] == "no_array"]
    lines += [f"- fragments (under {FRAGMENT_MAX_AA} aa): {', '.join(frag) or 'none'}",
              f"- not assessed (too little upstream sequence): {', '.join(na) or 'none'}",
              f"- assessed, no array: {', '.join(none) or 'none'}",
              "", "## Method",
              "Motif ORF search, then the paper's k-mer array scan (Methods p.32 step 1, "
              "100 mononucleotide shuffles) and array delimitation (step 2).",
              "", "## Limits",
              "Nine synthetic contigs and a motif in place of the 52 RT profile HMMs. "
              "No experiment here supports any biological claim."]
    return "\n".join(lines) + "\n"


def score(shared: Path, truth_path: Path) -> Scorecard:
    truth = json.loads(truth_path.read_text())
    cand = {r["contig"]: r["locus_id"] for r in _table(shared / "rt_candidates.tsv")}
    calls = {r["locus_id"]: r["status"] for r in _table(shared / "arrays.tsv")}
    sc = Scorecard()
    for t in truth:
        lid = cand.get(t["contig"])
        status = calls.get(lid) if lid else None
        if t["kind"] == "true_art":
            sc.true_arrays += 1
            if status == "array":
                sc.found_true += 1
            else:
                sc.missed.append(f"{t['contig']} ({status or 'never scanned'})")
        else:
            expected = {"no_array": "no_array", "short_upstream": "not_assessed",
                        "dense_decoy": "no_array", "fragment": None}[t["kind"]]
            sc.decoys_ok[f"{t['contig']} {t['kind']}"] = status == expected
            if status == "array":
                sc.false_calls.append(t["contig"])
    return sc


def _git_env() -> None:
    for k, v in (("GIT_AUTHOR_NAME", "art-demo"), ("GIT_COMMITTER_NAME", "art-demo"),
                 ("GIT_AUTHOR_EMAIL", "demo@localhost"), ("GIT_COMMITTER_EMAIL", "demo@localhost")):
        os.environ.setdefault(k, v)


def run_demo(opts: DemoOptions, say: Callable[[str], None] = print) -> DemoResult:
    out = opts.out.resolve()
    if out.exists():
        if not opts.force:
            raise FileExistsError(f"{out} exists; pass --force to replace it")
        shutil.rmtree(out)
    camp = out / "campaign"
    camp.mkdir(parents=True)
    narrator = Narrator(say, opts.explain)
    world_dir = out / "world"
    build_world(world_dir, opts.seed)
    genomes, truth = world_dir / "genomes.fna", world_dir / "truth.json"
    shared = camp / "shared"
    shared.mkdir()
    _git_env()
    subprocess.run(["git", "init", "-q", "-b", "campaign"], cwd=camp, check=False,
                   capture_output=True)

    if opts.agent == "ollama":
        check_ollama(opts.model, opts.ollama_url)
        agent = OllamaAgent(genomes, shared, OllamaChat(opts.model, opts.ollama_url))
        who = f"small model {opts.model} via Ollama"
    else:
        agent = RuleAgent(genomes, shared, opts.stage_revision)
        who = "rule policy (no model)"

    cfg = CampaignConfig(max_revisions=3, max_gate_failures=3, max_tasks_total=40,
                         dispatch_concurrency=opts.concurrency)
    store = RecordStore(camp)
    kb = KnowledgeBase(camp / "kb")
    ledger = SessionLedger(camp / "ledger.jsonl")
    backend = NarratedBackend(agent, narrator, store.records)
    orch = Orchestrator(cfg, store, kb, ledger, backend, make_gates(shared, genomes, narrator),
                        triage=make_triage(shared, narrator), skills_dir=None)

    say(f"ART harness demo | agents: {who} | world: 9 synthetic contigs, planted answer kept out of the tools")
    say("Five gated stages; each task = worker -> supervisor (-> revise) -> curator.\n")
    t0 = time.monotonic()
    error = None
    try:
        orch.run_stage_chain(BRIEFS)
    except RuntimeError as exc:
        error = str(exc)
        narrator.event("stop", why=error)
        say(f"\nCAMPAIGN STOPPED: {exc}")
    report = orch.report
    report.tasks_total = len(store.list_tasks())

    filed = None
    if error is None:
        synth = next((r for r in store.list_tasks()
                      if r.stage == STAGES[-1] and r.status is TaskStatus.CURATED), None)
        if synth:
            say("\nReport: draft from the shared tables, then the editor decides.")
            text = build_report(shared, [r.task_id for r in store.list_tasks()])
            filed = orch.file_report(synth.task_id, text)
            narrator.event("report", id=synth.task_id, filed=filed is not None)
            say(f"  report {'FILED at ' + str(filed.relative_to(out)) if filed else 'NOT filed'}",
                )
        report = orch.report
    sc = score(shared, truth) if (shared / "arrays.tsv").exists() else None
    wt = write_walkthrough(out, camp, store, ledger, report, sc, opts, who,
                           elapsed=time.monotonic() - t0, error=error, filed=filed)
    score_line = (f"{sc.found_true}/{sc.true_arrays} planted arrays found, "
                  f"{len(sc.false_calls)} false calls, decoys handled "
                  f"{sum(sc.decoys_ok.values())}/{len(sc.decoys_ok)}") if sc else ""
    write_replay(out, narrator.events, {
        "line": f"{who} \u00b7 {report.tasks_total} tasks \u00b7 {time.monotonic() - t0:.0f}s",
        "score": score_line, "summary": report.render().replace("\n", " | ")})
    say("\n" + report.render())
    if sc:
        say(f"Scorecard vs hidden truth: {sc.found_true}/{sc.true_arrays} planted arrays "
            f"found, {len(sc.false_calls)} false calls, decoys handled "
            f"{sum(sc.decoys_ok.values())}/{len(sc.decoys_ok)}")
    say(f"Walkthrough: {wt}\nReplay (open in a browser): {out / 'replay.html'}")
    return DemoResult(out, report, sc, error is None, error, filed, wt)
