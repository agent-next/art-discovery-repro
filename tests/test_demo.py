"""The demo campaign: a real end-to-end run of the workflow at toy scale.

Oracles: the world's planted truth (which the tools never read), the git history of
the task records, and the scripted checks. The small-model path is exercised against a
local HTTP stub that speaks the Ollama chat API, so nothing here needs a model.
"""

from __future__ import annotations

import json
import re
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest

from artharness.demo.agents import (
    OllamaAgent,
    OllamaChat,
    RuleAgent,
    check_ollama,
    evidence_problems,
    report_problems,
)
from artharness.demo.campaign import (
    DemoOptions,
    make_gates,
    make_triage,
    run_demo,
    score,
)
from artharness.demo.narrator import Narrator
from artharness.demo.tools import Sandbox, ToolError, run_tool
from artharness.demo.world import build_world
from artharness.runner.base import SessionSpec

# -- world and tools ----------------------------------------------------------------------


def test_world_is_deterministic_and_plants_the_advertised_answer(tmp_path: Path):
    a = build_world(tmp_path / "a", seed=1)
    build_world(tmp_path / "b", seed=1)
    assert (tmp_path / "a" / "genomes.fna").read_bytes() == \
        (tmp_path / "b" / "genomes.fna").read_bytes()
    kinds = sorted(t.kind for t in a)
    assert kinds.count("true_art") == 4
    assert {"short_upstream", "dense_decoy", "fragment", "no_array"} <= set(kinds)


@pytest.fixture(scope="module")
def sandbox(tmp_path_factory) -> Sandbox:
    root = tmp_path_factory.mktemp("sb")
    build_world(root / "w")
    (root / "art").mkdir()
    return Sandbox(root / "w" / "genomes.fna", root / "shared", root / "art")


def test_tools_recover_planted_loci_and_reject_decoys(sandbox: Sandbox):
    assert run_tool(sandbox, "list_contigs").headline.startswith("contigs=9 ")
    assert run_tool(sandbox, "find_rt_orfs").headline == "rt_orfs=9 contigs_with_rt=9"
    assert "fragments=1" in run_tool(sandbox, "classify_rt").headline
    scan = run_tool(sandbox, "scan_arrays")
    # the dense decoy (spacing < 100 nt) must NOT be called: the paper's suppression rule
    assert scan.headline == "assessed=8 arrays=4 no_array=3 not_assessed=1"
    assert len(scan.suggestions) == 4
    dive = run_tool(sandbox, "deep_dive", "locus=L01")
    assert "delimited=1" in dive.headline
    assert (sandbox.shared / "deep_dives" / "L01.tsv").exists()


def test_tools_refuse_bad_calls(sandbox: Sandbox, tmp_path: Path):
    with pytest.raises(ToolError, match="unknown tool"):
        run_tool(sandbox, "rm_rf")
    with pytest.raises(ToolError, match="no array call"):
        run_tool(sandbox, "deep_dive", "locus=L04")
    with pytest.raises(ToolError, match="bad arguments"):
        run_tool(sandbox, "find_rt_orfs", "scan_arrays")
    fresh = Sandbox(sandbox.genomes, tmp_path / "empty", tmp_path)
    with pytest.raises(ToolError, match="does not exist yet"):
        run_tool(fresh, "classify_rt")


# -- scripted checks (the supervisor's and editor's teeth) ---------------------------------

LOG = [{"tool": "scan_arrays", "args": "", "ok": True, "headline": "arrays=4 no_array=3",
        "output": "HEADLINE: arrays=4 no_array=3\nL01\tarray\tR=7"}]


def test_evidence_check_has_power():
    good = "# s\nHEADLINE: arrays=4 no_array=3\nFour calls."
    assert evidence_problems(good, LOG) == []
    assert any("no tool call" in p for p in evidence_problems(good, []))
    assert any("missing line" in p for p in evidence_problems("# s\nFour calls.", LOG))
    invented = evidence_problems(good + "\nCoverage was 87 percent.", LOG)
    assert any("87" in p for p in invented)
    assert evidence_problems(good + "\nLocus L01 has R=7.", LOG) == []


def test_editor_check_has_power(tmp_path: Path):
    shared = tmp_path
    (shared / "arrays.tsv").write_text("locus_id\tstatus\nL01\tarray\nL02\tno_array\n")
    (shared / "rt_candidates.tsv").write_text("locus_id\nL01\nL02\n")
    ok = "1 array-positive loci; L01 is the one, L02 is not."
    assert report_problems(ok, shared) == []
    assert any("L09" in p for p in report_problems(ok + " L09 too.", shared))
    assert any("2 array-positive" not in p and "must state" in p
               for p in report_problems("2 array-positive loci", shared))
    assert report_problems(ok, tmp_path / "nowhere")


def _tables(shared: Path, arrays: str) -> None:
    (shared / "contigs.tsv").write_text("contig\nc1\n")
    (shared / "rt_candidates.tsv").write_text(
        "locus_id\tcontig\tstrand\tstart\tend\taa_len\tmotif\nL01\tc1\t+\t10\t40\t10\tYADD\n"
        "L02\tc1\t+\t50\t80\t10\tYADD\n")
    (shared / "rt_classes.tsv").write_text(
        "locus_id\trt_class\tmotif\tfragment\nL01\tA\tYADD\t0\nL02\tA\tYADD\t0\n")
    (shared / "arrays.tsv").write_text(arrays)


def test_gates_fail_when_the_work_is_incomplete(tmp_path: Path):
    genomes = tmp_path / "g.fna"
    genomes.write_text(">c1\n" + "A" * 200 + "\n")
    shared = tmp_path / "shared"
    shared.mkdir()
    lines: list[str] = []
    gates = make_gates(shared, genomes, Narrator(lines.append, explain=False))
    assert not gates["1_input_assembly"]("x")  # nothing published yet
    _tables(shared, "locus_id\tstatus\nL01\tarray\n")  # L02 never scanned
    assert gates["1_input_assembly"]("x")
    assert gates["2_database_sweep"]("x")
    assert gates["3_rt_classification"]("x")
    assert not gates["4_neighborhood_census"]("x")
    assert "covers 1 of 2" in "\n".join(lines)
    _tables(shared, "locus_id\tstatus\nL01\tarray\nL02\tno_array\n")
    assert gates["4_neighborhood_census"]("x")
    assert not gates["5_deep_dives"]("x")  # L01 has an array but no deep dive
    (shared / "deep_dives").mkdir()
    (shared / "deep_dives" / "L01.tsv").write_text("x\n")
    assert gates["5_deep_dives"]("x")


def test_triage_rejects_with_a_reason(tmp_path: Path):
    shared = tmp_path
    (shared / "arrays.tsv").write_text(
        "locus_id\tstatus\nL01\tarray\nL02\tno_array\nL03\tnot_assessed\n")
    triage = make_triage(shared, Narrator(lambda _s: None, explain=False))

    class P:
        task_id = "t0001"

    assert triage("Deep dive locus L01", P()) == (True, "")
    ok, why = triage("Deep dive locus L01 again", P())
    assert not ok and "already queued" in why
    ok, why = triage("Deep dive locus L02", P())
    assert not ok and "'no_array'" in why
    assert "'not_assessed'" in triage("Deep dive locus L03", P())[1]
    assert "'absent'" in triage("Deep dive locus L07", P())[1]
    assert "names no locus" in triage("Validate expression levels", P())[1]


# -- the whole campaign, rule policy --------------------------------------------------------------


@pytest.fixture(scope="module")
def rule_run(tmp_path_factory):
    lines: list[str] = []
    res = run_demo(DemoOptions(out=tmp_path_factory.mktemp("demo") / "run"), lines.append)
    return res, "\n".join(lines)


def test_rule_campaign_runs_all_five_gated_stages_and_files_a_report(rule_run):
    res, console = rule_run
    assert res.error is None and res.filed is not None
    for stage in ("1_input_assembly", "2_database_sweep", "3_rt_classification",
                  "4_neighborhood_census", "5_deep_dives"):
        assert f"GATE {stage}: PASS" in console
    assert res.report.tasks_total == 9 and res.report.reports_filed == 1
    assert res.report.stalled == 0


def test_rule_campaign_is_scored_against_the_hidden_truth(rule_run):
    sc = rule_run[0].scorecard
    assert sc.perfect and sc.found_true == sc.true_arrays == 4


def test_the_staged_sloppy_pass_is_revised_then_accepted(rule_run):
    res = rule_run[0]
    t4 = res.out / "campaign" / "records" / "t0004"
    assert "no tool call succeeded" in (t4 / "verdict-r1.md").read_text()
    assert (t4 / "verdict.md").read_text().startswith("VERDICT: accept")
    assert res.report.revised == 1


def test_the_supervisors_overeager_followup_is_rejected_with_a_reason(rule_run):
    res, console = rule_run
    assert res.report.rejected_at_triage == 1
    rej = next((res.out / "campaign" / "records" / "t0004").glob("triage-rejection-*.md"))
    assert "L06" in rej.read_text() and "not_assessed" in rej.read_text()


def test_task_state_changes_are_git_commits_in_order(rule_run):
    camp = rule_run[0].out / "campaign"
    log = subprocess.run(["git", "log", "--reverse", "--format=%s"], cwd=camp,
                         capture_output=True, text=True, check=True).stdout.splitlines()
    idx = {m: i for i, m in enumerate(log)}
    assert idx["task(t0004): supervisor revise (1/3)"] < idx["task(t0004): supervisor accepted"]
    assert idx["task(t0001): curated"] < idx["task(t0002): open  [seed]"]  # gate order
    assert any(m.startswith("report(t0009): filed") for m in log)


def test_walkthrough_is_generated_from_the_run(rule_run):
    res = rule_run[0]
    text = res.walkthrough.read_text()
    for needle in ("## Why the workflow is shaped this way", "## A revision, in full",
                   "## Triage rejections", "planted arrays found: **4/4**",
                   "## Audit trail", "supervisor revise (1/3)"):
        assert needle in text


def _replay_data(run_dir: Path) -> dict:
    page = (run_dir / "replay.html").read_text()
    m = re.search(r'<script id="data" type="application/json">(.*?)</script>', page, re.S)
    assert m, "replay page lost its embedded data"
    return json.loads(m.group(1))


def test_replay_events_agree_with_the_campaign_report(rule_run):
    res = rule_run[0]
    events = _replay_data(res.out)["events"]
    kinds = [e["k"] for e in events]
    assert json.loads((res.out / "events.json").read_text()) == events
    assert kinds.count("task") == res.report.tasks_total
    assert kinds.count("gate") == 5 and all(e["ok"] for e in events if e["k"] == "gate")
    assert sum(1 for e in events if e["k"] == "triage" and not e["ok"]) == 1
    assert sum(1 for e in events if e["k"] == "session" and e["verdict"] == "revise") == 1
    assert kinds[-1] == "report" and events[-1]["filed"] is True


def test_replay_page_is_self_contained_and_survives_hostile_text(tmp_path: Path):
    from artharness.demo.replay import write_replay

    nasty = "</script><script>alert(1)</script>"
    path = write_replay(tmp_path, [{"k": "task", "id": "t1", "brief": nasty}],
                        {"line": nasty, "score": "", "summary": ""})
    page = path.read_text()
    assert page.count("</script>") == 2  # the data block and the player, nothing injected
    assert not re.search(r'(src|href)="', page) and "https://" not in page
    assert _replay_data(tmp_path)["events"][0]["brief"] == nasty


def test_replay_stations_do_not_overlap_and_their_text_fits():
    from artharness.demo.replay import TEMPLATE

    js = TEMPLATE.read_text()
    block = re.search(r"const ST = \{(.*?)\};", js, re.S).group(1)
    st = {k: tuple(map(int, v.split(","))) for k, v in
          re.findall(r"(\w+):\[([\d,]+)\]", block)}
    assert len(st) == 10
    view_w, view_h = map(int, re.search(r'viewBox="0 0 (\d+) (\d+)"', js).groups())
    for k, (x, y, w, h) in st.items():
        assert x >= 0 and y >= 0 and x + w <= view_w and y + h <= view_h, k
    names = sorted(st)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ax, ay, aw, ah = st[a]
            bx, by, bw, bh = st[b]
            assert not (ax < bx + bw and bx < ax + aw and ay < by + bh and by < ay + ah), (a, b)
    # 18px bold titles are ~10.5 px per character; sub-lines wrap at 21 characters of 13.5px text
    titles = dict(re.findall(r"(\w+):\['([^']+)'", re.search(r"const TITLE = \{(.*?)\};", js, re.S).group(1)))
    for k, title in titles.items():
        assert len(title) * 10.5 <= st[k][2] - 8, (k, title)


def test_readme_hero_svg_is_well_formed_and_its_timeline_is_wired():
    import xml.etree.ElementTree as ET

    svg = Path(__file__).parent.parent / "docs" / "assets" / "how-it-works.svg"
    root = ET.parse(svg).getroot()
    ids = {el.get("id") for el in root.iter() if el.get("id")}
    refs = set(re.findall(r'begin="[^"]*?\b(s\d+)\.(?:begin|end)', svg.read_text()))
    assert {f"s{i}" for i in range(1, 11)} <= ids and refs <= ids
    assert not any(el.tag.endswith("script") for el in root.iter())


def test_no_stage_revision_flag_removes_the_revision(tmp_path: Path):
    res = run_demo(DemoOptions(out=tmp_path / "r", stage_revision=False, explain=False),
                   lambda _s: None)
    assert res.report.revised == 0 and res.error is None


def test_existing_output_dir_needs_force(rule_run):
    with pytest.raises(FileExistsError):
        run_demo(DemoOptions(out=rule_run[0].out), lambda _s: None)


def test_score_marks_a_wrong_call(rule_run, tmp_path: Path):
    res = rule_run[0]
    shared = tmp_path / "shared"
    shared.mkdir()
    for f in (res.out / "campaign" / "shared").glob("*.tsv"):
        (shared / f.name).write_text(f.read_text())
    text = (shared / "arrays.tsv").read_text().replace("L04\tno_array", "L04\tarray")
    (shared / "arrays.tsv").write_text(text)
    sc = score(shared, res.out / "world" / "truth.json")
    assert not sc.perfect and sc.false_calls == ["phage_04"]


# -- the small-model path against an Ollama stub ---------------------------------------------------


class Stub:
    """Speaks /api/chat and /api/tags. ``respond(messages, options)`` returns text."""

    def __init__(self, respond, models=("qwen3:0.6b",)):
        self.respond, self.models, self.requests = respond, models, []
        outer = self

        class H(BaseHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def do_GET(self):
                body = json.dumps({"models": [{"name": m} for m in outer.models]}).encode()
                self.send_response(200)
                self.end_headers()
                self.wfile.write(body)

            def do_POST(self):
                req = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
                outer.requests.append(req)
                text = outer.respond(req["messages"], req.get("options", {}))
                body = json.dumps({"message": {"content": text}, "prompt_eval_count": 11,
                                   "eval_count": 5}).encode()
                self.send_response(200)
                self.end_headers()
                self.wfile.write(body)

        self.srv = HTTPServer(("127.0.0.1", 0), H)
        threading.Thread(target=self.srv.serve_forever, daemon=True).start()
        self.url = f"http://127.0.0.1:{self.srv.server_port}"

    def close(self):
        self.srv.shutdown()


def compliant(messages, options):
    """A well-behaved model: follows the protocol it is shown."""
    last = messages[-1]["content"]
    system = messages[0]["content"]
    if "RESULT of" in last:
        heads = re.findall(r"^HEADLINE: .*$", last, re.MULTILINE)
        sugg = re.findall(r"^SUGGEST_FOLLOWUP: (.*)$", last, re.MULTILINE)
        return "SUMMARY:\n" + "\n".join(heads) + "\nDone.\n" + "\n".join(
            f"PROPOSE_FOLLOWUP: {s}" for s in sugg)
    if "You are the supervisor" in system:
        return "Looks right.\nVERDICT: accept"
    if "You are the curator" in system:
        return "Findings recorded."
    if "You are the editor" in system:
        return "FILE: yes"
    call = re.search(r"line `TOOL: ([^`]+)`", last)
    return "PLAN: do it\n- TOOL: " + (call.group(1) if call else "list_contigs")


def _spec(role: str, workdir: Path, prompt: str, brief: str) -> SessionSpec:
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "brief.md").write_text(brief)
    system = {"worker": "You are the worker", "supervisor": "You are the supervisor",
              "curator": "You are the curator", "editor": "You are the editor"}[role]
    return SessionSpec(role=role, task_id="t0001", system_prompt=system, user_prompt=prompt,
                       workdir=workdir)


@pytest.fixture()
def agent_env(tmp_path: Path):
    build_world(tmp_path / "w")
    shared = tmp_path / "shared"
    shared.mkdir()
    made = []

    def make(respond):
        stub = Stub(respond)
        made.append(stub)
        return stub, OllamaAgent(tmp_path / "w" / "genomes.fna", shared,
                                 OllamaChat("qwen3:0.6b", stub.url))

    yield make, tmp_path
    for s in made:
        s.close()


BRIEF = "Stage 1 - inventory.\nTools: list_contigs"
PROMPT = "# Task t0001 (1_input_assembly)\n\n" + BRIEF


def test_ollama_worker_runs_tools_and_writes_a_checkable_summary(agent_env):
    make, tmp = agent_env
    stub, agent = make(compliant)
    wd = tmp / "rec" / "t0001"
    out = agent.run(_spec("worker", wd, PROMPT, BRIEF))
    assert "HEADLINE: contigs=9" in (wd / "summary.md").read_text()
    assert agent.evidence_review(_spec("supervisor", wd, PROMPT, BRIEF)) == []
    assert out.result.input_tokens_uncached == 22 and out.result.output_tokens == 10
    assert "TOOL: list_contigs" in stub.requests[0]["messages"][-1]["content"]


def test_ollama_worker_cannot_leave_its_task_scope(agent_env):
    make, tmp = agent_env
    stub, agent = make(lambda m, o: "PLAN: x\nTOOL: list_contigs\nTOOL: deep_dive locus=L01"
                       if len(m) == 2 else "SUMMARY:\nHEADLINE: contigs=9 total_nt=1")
    wd = tmp / "rec" / "t0001"
    agent.run(_spec("worker", wd, PROMPT, BRIEF))
    log = [json.loads(ln) for ln in (wd / "artifacts" / "tool_log.jsonl").read_text().splitlines()]
    assert [(e["tool"], e["ok"]) for e in log] == [("list_contigs", True), ("deep_dive", False)]
    assert "outside this task's scope" in log[1]["output"]


def test_ollama_worker_invented_results_are_cut_and_run_nothing(agent_env):
    make, tmp = agent_env
    stub, agent = make(lambda m, o: "PLAN: x\nRESULT: 9 contigs, all fine\nSUMMARY: done")
    wd = tmp / "rec" / "t0001"
    agent.run(_spec("worker", wd, PROMPT, BRIEF))
    assert not (wd / "artifacts" / "tool_log.jsonl").exists()
    sup = agent.evidence_review(_spec("supervisor", wd, PROMPT, BRIEF))
    assert any("no tool call succeeded" in p for p in sup)


def test_ollama_revision_reuses_tool_results_and_moves_the_seed(agent_env):
    make, tmp = agent_env
    stub, agent = make(compliant)
    wd = tmp / "rec" / "t0001"
    agent.run(_spec("worker", wd, PROMPT, BRIEF))
    first = json.loads((wd / "artifacts" / "tool_log.jsonl").read_text().splitlines()[0])
    (wd / "verdict-r1.md").write_text("VERDICT: revise\n- quote the headline\n")
    n_before = len(stub.requests)
    agent.run(_spec("worker", wd, PROMPT, BRIEF))
    assert len(stub.requests) == n_before + 1  # one turn: the summary, no new tool round
    lines = (wd / "artifacts" / "tool_log.jsonl").read_text().splitlines()
    assert len(lines) == 1 and json.loads(lines[0]) == first
    assert stub.requests[-1]["options"]["seed"] == stub.requests[0]["options"]["seed"] + 1


def test_ollama_supervisor_verdict_needs_both_the_script_and_the_model(agent_env):
    make, tmp = agent_env
    wd = tmp / "rec" / "t0001"
    stub, agent = make(compliant)
    agent.run(_spec("worker", wd, PROMPT, BRIEF))
    sup = _spec("supervisor", wd, PROMPT, BRIEF)
    assert agent.run(sup).verdict == "accept"
    # model says revise although the script passes: honoured
    _, strict = make(lambda m, o: "VERDICT: revise")
    assert strict.run(sup).verdict == "revise"
    # model rubber-stamps although the script fails: overridden
    (wd / "summary.md").write_text("# s\nAll good, 99 percent.\n")
    _, lenient = make(lambda m, o: "Great work.\nVERDICT: accept")
    out = lenient.run(sup)
    assert out.verdict == "revise" and "missing line" in out.verdict_notes
    # no parseable verdict: None, which the orchestrator treats as revise
    _, mute = make(lambda m, o: "hmm")
    (wd / "summary.md").write_text(f"# s\nHEADLINE: contigs=9 total_nt={_total(wd)}\n")
    assert mute.run(sup).verdict is None


def _total(wd: Path) -> str:
    return re.search(r"total_nt=(\d+)", (wd / "artifacts" / "tool_log.jsonl").read_text()).group(1)


def test_ollama_editor_cannot_file_a_report_the_script_rejects(agent_env):
    make, tmp = agent_env
    stub, agent = make(lambda m, o: "FILE: yes")
    spec = _spec("editor", tmp / "ed", "9 array-positive loci L01", "")
    out = agent.run(spec)  # no tables at all
    assert out.verdict == "no"


def test_ollama_curator_returns_prose(agent_env):
    make, tmp = agent_env
    stub, agent = make(lambda m, o: "Nine contigs were listed.")
    assert agent.run(_spec("curator", tmp / "c", PROMPT, BRIEF)).text == "Nine contigs were listed."


def test_full_campaign_with_a_compliant_stub_model(tmp_path: Path):
    stub = Stub(compliant)
    try:
        lines: list[str] = []
        res = run_demo(DemoOptions(out=tmp_path / "r", agent="ollama", ollama_url=stub.url),
                       lines.append)
    finally:
        stub.close()
    assert res.error is None and res.filed is not None
    assert res.scorecard.perfect
    assert res.report.stalled == 0 and res.report.reports_filed == 1
    assert "small model qwen3:0.6b via Ollama" in "\n".join(lines)
    assert res.report.follow_ups == 4  # the model copied the SUGGEST lines into follow-ups


def test_a_model_that_never_uses_tools_is_stopped_by_the_gates(tmp_path: Path):
    stub = Stub(lambda m, o: "I have thought about it.")
    try:
        res = run_demo(DemoOptions(out=tmp_path / "r", agent="ollama", ollama_url=stub.url,
                                   explain=False), lambda _s: None)
    finally:
        stub.close()
    assert res.error and "gate for 1_input_assembly failed" in res.error
    assert res.report.tasks_total == 1 and res.report.stalled == 1
    assert "stopped early" in res.walkthrough.read_text()


def test_model_server_checks(tmp_path: Path, monkeypatch):
    stub = Stub(compliant)
    try:
        check_ollama("qwen3:0.6b", stub.url)
        with pytest.raises(RuntimeError, match="not installed"):
            check_ollama("llama9", stub.url)
    finally:
        stub.close()
    with pytest.raises(RuntimeError, match="cannot reach"):
        check_ollama("qwen3:0.6b", "http://127.0.0.1:9")
    monkeypatch.delenv("ARTHARNESS_ALLOW_LIVE", raising=False)
    with pytest.raises(RuntimeError, match="ARTHARNESS_ALLOW_LIVE"):
        check_ollama("m", "http://models.example.com:11434")


def test_rule_agent_reads_its_brief_not_a_hardcoded_script(tmp_path: Path):
    build_world(tmp_path / "w")
    shared = tmp_path / "shared"
    shared.mkdir()
    agent = RuleAgent(tmp_path / "w" / "genomes.fna", shared, stage_revision=False)
    wd = tmp_path / "rec" / "t0001"
    agent.run(_spec("worker", wd, PROMPT, BRIEF))
    assert (shared / "contigs.tsv").exists() and not (shared / "rt_candidates.tsv").exists()
