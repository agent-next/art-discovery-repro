"""WALKTHROUGH.md: the finished run explained from its own records.

Everything in the file is read back from the campaign directory (task records, git
history, knowledge base, ledger, tables), so it describes what this run did, not what
the harness is supposed to do.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from ..accounting import SessionLedger
from ..orchestrator import STAGES, CampaignReport
from ..records import RecordStore

PAPER = {"tasks": 119, "sessions": 949, "revised": 49, "rejected": 10, "reports": 19}

STAGE_WHY = {
    "1_input_assembly": "know exactly what is being analysed before anything else",
    "2_database_sweep": "cast the widest net first (paper: 52 RT profile HMMs over 1.9B clusters)",
    "3_rt_classification": "sort candidates so later stages compare like with like "
                           "(paper: 9 classes)",
    "4_neighborhood_census": "look next to every RT for an array (paper: 10,983 loci)",
    "5_deep_dives": "spend effort only on what the census flagged (paper: 16 deep dives)",
}


def _log(camp: Path) -> list[str]:
    try:
        r = subprocess.run(["git", "log", "--reverse", "--format=%s"], cwd=camp,
                           capture_output=True, text=True, check=True)
        return r.stdout.splitlines()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []


def _headlines(store: RecordStore, tid: str) -> str:
    p = store.artifact_dir(tid) / "tool_log.jsonl"
    if not p.exists():
        return "-"
    heads = [json.loads(ln)["headline"] for ln in p.read_text().splitlines() if ln.strip()]
    return "; ".join(h for h in heads if h) or "(tool errors only)"


def write_walkthrough(out: Path, camp: Path, store: RecordStore, ledger: SessionLedger,
                      report: CampaignReport, sc, opts, who: str, elapsed: float,
                      error: str | None, filed: Path | None) -> Path:
    tasks = store.list_tasks()
    summ = ledger.summary()
    L: list[str] = []
    w = L.append
    w("# ART harness demo: what happened, and why")
    w("")
    w(f"Agents: **{who}**. Wall time {elapsed:.0f}s. Everything below was read back from "
      f"`{camp.name}/` after the run.")
    if error:
        w(f"\n**The campaign stopped early: {error}.** That is the harness working as "
          "designed: a failed scripted gate keeps later stages closed.")
    w("")
    w("## What this is")
    w("A miniature of the campaign in Yoon et al. 2026 (Anthropic). The paper had "
      "~950 model sessions search 1.9B protein clusters and find ART, a tandem-repeat-"
      "array reverse transcriptase system. This run keeps the paper's **workflow** and "
      "shrinks the **world** to nine synthetic contigs with a hidden, planted answer, "
      "and the **model** to a rule policy or a 0.6B local model.")
    w("")
    w("## Why the workflow is shaped this way")
    w("| Harness element | What it does | Why (paper) | Seen in this run |")
    w("| --- | --- | --- | --- |")
    w("| Stage chain + scripted gates | five stages; a stage's tasks open only after the "
      "previous gate passes | model output is not trusted to say \"done\" (p.28) | "
      f"{len(STAGES)} gates, all checked by code |")
    w("| Worker | plans, runs tools, writes summary | the agent that does the science | "
      f"{summ['sessions_by_role']['worker']} sessions |")
    w("| Supervisor | accepts or returns for revision after checking claims against "
      "artifacts | independent review; 49 of 119 paper tasks were revised | "
      f"{summ['sessions_by_role']['supervisor']} sessions, {report.revised} task(s) revised |")
    w("| Curator + knowledge base | writes findings into a store later agents read | "
      "agents do not repeat or contradict earlier work | "
      f"{summ['sessions_by_role']['curator']} entries |")
    w("| Triage | follow-ups are released or rejected with a written reason | keeps "
      "agent-proposed work bounded | "
      f"{report.follow_ups} released, {report.rejected_at_triage} rejected |")
    w("| Editor | reviews the report before it is filed | unsupported claims are caught "
      "before they leave | " + ("report filed" if filed else "report not filed") + " |")
    w("| Git task records | every state change is a commit | the audit trail humans "
      "monitor | see the log below |")
    w("")
    w("## The five stages")
    w("| stage | why it exists | tasks | status |")
    w("| --- | --- | --- | --- |")
    for s in STAGES:
        ts = [t for t in tasks if t.stage == s]
        w(f"| {s} | {STAGE_WHY[s]} | {', '.join(t.task_id for t in ts) or '-'} | "
          f"{', '.join(sorted({t.status.value for t in ts})) or 'not opened'} |")
    w("")
    w("## Every task")
    w("| task | stage | origin | revisions | status | tool headline |")
    w("| --- | --- | --- | --- | --- | --- |")
    for t in tasks:
        origin = t.origin.value + (f" of {t.parent}" if t.parent else "")
        w(f"| {t.task_id} | {t.stage} | {origin} | {t.revisions} | {t.status.value} | "
          f"`{_headlines(store, t.task_id)}` |")
    w("")
    revised = [t for t in tasks if t.revisions]
    if revised:
        w("## A revision, in full")
        t = revised[0]
        w(f"Task {t.task_id} was returned {t.revisions} time(s). The supervisor's first "
          f"reason (`{t.task_id}/verdict-r1.md`):")
        w("")
        w("```")
        w((store.records / t.task_id / "verdict-r1.md").read_text().strip())
        w("```")
        w("")
    rej = []
    for p in sorted(store.records.glob("*/triage-rejection-*.md")):
        rej.append(p.read_text().strip())
    lp = camp / "triage-rejections.log"
    if lp.exists():
        rej.append(lp.read_text().strip())
    if rej:
        w("## Triage rejections (each carries a written reason)")
        for r in rej:
            w("```")
            w(r)
            w("```")
        w("")
    w("## Knowledge base")
    for path, text in KB_ENTRIES(camp):
        body = text.split("---", 2)[-1].strip().replace("\n", " ")
        w(f"- `{path}`: {body[:220]}")
    w("")
    w("## Accounting next to the paper")
    w("| | this run | paper |")
    w("| --- | --- | --- |")
    w(f"| tasks | {report.tasks_total} | {PAPER['tasks']} |")
    w(f"| sessions | {summ['sessions_total']} | {PAPER['sessions']} |")
    w(f"| tasks revised >= 1x | {report.revised} | {PAPER['revised']} |")
    w(f"| triage rejections | {report.rejected_at_triage} | {PAPER['rejected']} |")
    w(f"| reports filed | {report.reports_filed} | {PAPER['reports']} |")
    tok = summ["tokens"]
    w(f"| model tokens (in/out) | {tok['input_uncached']}/{tok['output']} | "
      "11.3M/14.9M (+189.5M cache writes) |")
    w("")
    if sc:
        w("## Scored against the hidden truth")
        w(f"- planted arrays found: **{sc.found_true}/{sc.true_arrays}**")
        w(f"- false array calls: **{len(sc.false_calls)}** {sc.false_calls or ''}")
        if sc.missed:
            w(f"- missed: {sc.missed}")
        for k, v in sc.decoys_ok.items():
            w(f"- decoy `{k}`: {'handled correctly' if v else 'MISHANDLED'}")
        w("")
    if filed:
        w("## The filed report")
        w("```markdown")
        w(filed.read_text().strip())
        w("```")
        w("")
    w("## Audit trail (git log of the task records, oldest first)")
    w("```")
    L.extend(_log(camp))
    w("```")
    w("")
    w("## What is real and what is a stand-in")
    w("| real | stand-in |")
    w("| --- | --- |")
    w("| the orchestrator, gates, triage, roles, records, ledger | 9 synthetic contigs |")
    w("| the paper's k-mer array scan and array delimitation (`artharness.arrays`) | an "
      "RT motif instead of the 52 profile HMMs |")
    w("| the supervisor's claim check (code) | " + who + " instead of Claude Mythos 5 |")
    w("")
    w("## Try it")
    w("```")
    w("python -m artharness.demo                      # rule policy, offline")
    w("python -m artharness.demo --agent ollama      # small local model (qwen3:0.6b)")
    w("python -m artharness.demo --no-stage-revision # skip the staged sloppy pass")
    w("```")
    path = out / "WALKTHROUGH.md"
    path.write_text("\n".join(L) + "\n")
    return path


def KB_ENTRIES(camp: Path) -> list[tuple[str, str]]:
    d = camp / "kb" / "entries"
    return [(p.name, p.read_text()) for p in sorted(d.glob("*.md"))]
