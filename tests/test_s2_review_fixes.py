"""S2 regression-power tests for the /my-review S1 finder round (2026-09-24).

Every test here corresponds to a confirmed S1 finding; each must fail on the
pre-fix code (disabled-vs-enabled power per REVIEW.md bars).
"""

import random
from pathlib import Path

from artharness.arrays import (
    _chains,
    _regular_runs,
    build_pwm,
    delimit_array,
    pwm_extend,
)
from artharness.config import CampaignConfig
from artharness.records import RecordStore, TaskOrigin, TaskStatus
from artharness.runner.base import ScriptedBackend

# ---------------------------------------------------------------------------
# B1: delimitation skipped-copy rule anchored on the running median
# ---------------------------------------------------------------------------

def test_chains_allow_a_single_skipped_copy_mid_run():
    # copies at 0,200,400,800,1000: the 400->800 gap (2x median) must be
    # tolerated as one skipped copy -> a 5-copy chain; the old outside-range
    # test made the skip branch unreachable and broke the chain at 800.
    chain = _chains([0, 200, 400, 800, 1000])
    assert any(len(c) == 5 for c in chain), chain


def test_chains_reject_a_second_skipped_copy():
    # gaps 200,200,400(one skip -> legal 4-copy chain),400: growing to 5 would
    # consume a SECOND skipped copy, which the paper does not allow.
    chain = _chains([0, 200, 400, 800, 1200])
    assert any(len(c) == 4 for c in chain), chain  # single skip applies
    assert not any(len(c) == 5 for c in chain), chain  # second skip refused


# ---------------------------------------------------------------------------
# B2: one irregular gap must not discard an otherwise regular run
# ---------------------------------------------------------------------------

def test_regular_runs_keep_the_longest_tolerant_window():
    runs = _regular_runs([0, 200, 400, 700, 900, 1100])
    # the whole run fails the median check (gap 300 vs median 200), but the
    # two 3-copy halves are each regular; the old code returned [].
    assert runs and len(runs[0]) == 3


# ---------------------------------------------------------------------------
# B4: true median -- a legal 3-copy array must pass the tolerance rule
# ---------------------------------------------------------------------------

def test_regular_runs_accept_three_copies_with_gap_ratio_above_one_half():
    # gaps 200/350 (ratio 0.57): both within 30% of the true median 275.
    runs = _regular_runs([0, 200, 550])
    assert runs and len(runs[0]) == 3


# ---------------------------------------------------------------------------
# B5/B6: PWM edge cases
# ---------------------------------------------------------------------------

def test_pwm_extend_never_emits_negative_copy_starts():
    rng = random.Random(0)
    word = "ACGTACGTTGCAAT"
    # a PWM block hit in the first `off` positions of the window would have
    # produced copy starts before the window itself
    from artharness.arrays import DelimitedArray
    upstream = word * 60
    arr = DelimitedArray(locus="t", copy_starts=[100, 350, 600],
                         repeat=word[:10], score=5.0, shuffles_used=200,
                         spacings=[250, 250], block_offset=3)
    out = pwm_extend(arr, upstream, rng)
    assert all(p >= 0 for p in out.copy_starts)


def test_build_pwm_survives_a_zero_frequency_background_base(tmp_path: Path):
    # window with no T at all: the old dict lookup divided by background 0.0
    background = {"A": 0.5, "C": 0.5, "G": 0.0, "T": 0.0}
    pwm = build_pwm(["ACGTACGT", "ACGTACGT"], background)
    assert len(pwm) == 8


# ---------------------------------------------------------------------------
# B-minor-9 + delimitation sanity: delimit_array still returns on planted input
# ---------------------------------------------------------------------------

def test_delimit_array_still_recovers_a_planted_array():
    rng = random.Random(7)
    copies = []
    for _ in range(5):
        copies.append("ACGTACGTTA")
        copies.append("".join(rng.choice("AT") for _ in range(50)))
    upstream = "".join(copies)  # 300 nt, copies every 60 nt
    arr = delimit_array("locus", upstream, rng)
    assert arr is not None and len(arr.copy_starts) >= 4


# ---------------------------------------------------------------------------
# A1/A3/D2/D-4: backend parsing
# ---------------------------------------------------------------------------

def _backend_output(text, role="supervisor"):
    from artharness.runner.base import ClaudeCodeBackend, SessionSpec

    class FakeConfig:
        live_backend_env = "ARTHARNESS_ALLOW_LIVE"

    import os
    os.environ["ARTHARNESS_ALLOW_LIVE"] = "1"
    try:
        backend = ClaudeCodeBackend("test-model", FakeConfig())
    finally:
        del os.environ["ARTHARNESS_ALLOW_LIVE"]

    class FakeProc:
        stdout = "{}"

    spec = SessionSpec(role=role, task_id="t0001", system_prompt="s",
                       user_prompt="u", workdir=Path("."))
    # parse() is exercised through a stubbed subprocess call
    import json as _json

    payload = _json.dumps({"result": text, "usage": {"input_tokens": 1,
                                                     "output_tokens": 2}})

    class P:
        stdout = payload

    orig = __import__("subprocess").run
    __import__("subprocess").run = lambda *a, **k: P()
    try:
        out = backend.run(spec)
    finally:
        __import__("subprocess").run = orig
    return out


def test_verdict_uses_the_last_match_not_the_first():
    # a restated format example in the middle of the text must not beat the
    # actual verdict at the end (A3)
    text = "VERDICT: accept\n\n(format: VERDICT: revise)\n\nVERDICT: revise\nnotes"
    assert _backend_output(text).verdict == "revise"


def test_editor_file_uses_the_last_match():
    text = "FILE: no\n(final answer)\nFILE: yes"
    assert _backend_output(text, role="editor").verdict == "yes"


def test_supervisor_task_briefs_are_parsed_as_blocks():
    text = ("VERDICT: accept\n\nPROPOSED_TASK_BRIEF: scan Madawaska upstream "
            "region for the partner ORF.\n\nVERDICT notes end.")
    out = _backend_output(text)
    assert out.proposed_followups, "supervisor brief must be parsed (A1)"
    assert "Madawaska" in out.proposed_followups[0]


# ---------------------------------------------------------------------------
# A7: revision passes pass the completion check (stale summary withdrawn)
# ---------------------------------------------------------------------------

class EmptyRevision(ScriptedBackend):
    """Pass 1 writes summary.md; later worker passes produce nothing (the
    scripted writer is bypassed by deleting AFTER it ran — S3 finding 5:
    deleting before let it resurrect the file and the check was vacuous)."""

    def __init__(self):
        super().__init__()
        self.worker_calls = 0

    def run(self, spec):
        out = ScriptedBackend.run(self, spec)
        if spec.role == "worker":
            self.worker_calls += 1
            if self.worker_calls >= 2:
                (spec.workdir / "summary.md").unlink(missing_ok=True)
        if spec.role == "supervisor":
            out.verdict = "revise"  # never accepts; stall must come from checks
        return out


def test_revision_pass_failing_completion_check_stalls(tmp_path: Path):
    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import Orchestrator

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    kb = KnowledgeBase(root / "kb")
    ledger = SessionLedger(root / "ledger.jsonl")
    cfg = CampaignConfig(max_revisions=10, max_gate_failures=2)
    orch = Orchestrator(cfg, store, kb, ledger, EmptyRevision(),
                        {s: (lambda s: True) for s in
                         ("1_input_assembly", "2_database_sweep",
                          "3_rt_classification", "4_neighborhood_census",
                          "5_deep_dives")})
    orch.run_stage_chain({"1_input_assembly": ["do a thing"]})
    rec = store.list_tasks()[0]
    final = store.get(rec.task_id)
    # the revision pass produced no summary -> completion-check failure -> STALL.
    # Assert the COUNTERS, not just the status: a revisions-limit stall also
    # ends STALLED, which is exactly how the vacuous version passed (S3 no. 5).
    assert final.status == TaskStatus.STALLED
    # counters pin the MODE: the supervisor revised once, then both empty
    # revision passes failed the completion check (the pre-S3-fix code reached
    # revisions=10 / gate_failures=9 here — the check path was unreachable)
    assert final.gate_failures == 2, final.gate_failures
    assert final.revisions == 1, final.revisions
    assert orch.report.stalled == 1


# ---------------------------------------------------------------------------
# A5: triage rejection filenames do not collide
# ---------------------------------------------------------------------------

def test_two_rejections_same_parent_get_distinct_files(tmp_path: Path):
    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import Orchestrator

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    kb = KnowledgeBase(root / "kb")
    ledger = SessionLedger(root / "ledger.jsonl")
    backend = ScriptedBackend()
    orch = Orchestrator(CampaignConfig(), store, kb, ledger, backend,
                        {s: (lambda s: True) for s in
                         ("1_input_assembly", "2_database_sweep",
                          "3_rt_classification", "4_neighborhood_census",
                          "5_deep_dives")},
                        triage=lambda brief, parent: (False, "out of scope"))
    orch.run_stage_chain({"1_input_assembly": ["task one"]})
    parent = store.list_tasks()[0]
    orch.propose_followup("follow up A", parent)
    orch.propose_followup("follow up B", parent)
    files = sorted(p.name for p in
                   (store.records / parent.task_id).glob("triage-rejection-*.md"))
    assert len(files) == 2, files  # the old naming overwrote the first reason


# ---------------------------------------------------------------------------
# A4/A9: budget-capped seeding + tasks_total accounting
# ---------------------------------------------------------------------------

def test_stage_chain_survives_the_task_budget(tmp_path: Path):
    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import STAGES, Orchestrator

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    kb = KnowledgeBase(root / "kb")
    ledger = SessionLedger(root / "ledger.jsonl")
    orch = Orchestrator(CampaignConfig(max_tasks_total=2), store, kb, ledger,
                        ScriptedBackend(),
                        {s: (lambda s: True) for s in STAGES})
    orch.run_stage_chain({s: [f"brief {s}", f"brief {s} 2"] for s in STAGES})
    assert len(store.list_tasks()) == 2  # old code raised AttributeError here


# ---------------------------------------------------------------------------
# D-1: the soundness auto-lose rule must have regression power
# ---------------------------------------------------------------------------

def test_soundness_auto_lose_is_load_bearing():

    from artharness.config import CampaignConfig
    from artharness.tournament import JudgeScores, run_tournament

    # B's weighted total BEATS A's; only the soundness<=2 rule stops B winning.
    strong_unsound = JudgeScores(impact=5, novelty=5, soundness=2, actionability=5)
    reports = {"a": "A", "b": "B"}
    judge = lambda name, a, other, b: (  # noqa: E731
        strong_unsound if a == "B" else JudgeScores(4, 4, 5, 3))

    out = run_tournament(reports, judge, cfg=CampaignConfig(),
                         rng=random.Random(0))
    assert out["ranking"][0] == "a"  # auto-lose flips the outcome
    # power: with the rule disabled, B must win on weighted score
    wa = 4 * 0.35 + 4 * 0.30 + 5 * 0.25 + 3 * 0.10
    wb = 5 * 0.35 + 5 * 0.30 + 2 * 0.25 + 5 * 0.10
    assert wb > wa  # 4.25 vs 4.15: the rule, not the score, decides


def test_exact_ties_and_double_auto_lose_are_discarded_games():
    # Audit 2026-09-28 finding 9: an exact weighted tie (and a double
    # soundness auto-lose) recorded NO winner, but bradley_terry hard-coded
    # n_ij=2 — the fit assumed a decisive game that never happened. Ties must
    # be discarded games with per-pair n_ij accounting (INTERPRETED: the
    # paper states only "the report with the higher weighted score won").
    from artharness.config import CampaignConfig
    from artharness.tournament import JudgeScores, run_tournament

    scores = JudgeScores(4, 4, 4, 4)  # every game an exact tie
    reports = {"a": "A", "b": "B", "c": "C"}
    judge = lambda *args: scores  # noqa: E731

    out = run_tournament(reports, judge, cfg=CampaignConfig(),
                         rng=random.Random(0))
    assert out["games"] == 6  # every ordered pair was judged
    assert out["discarded_games"] == 6  # none produced a winner
    assert sum(out["wins"].values()) == 0
    # zero decisive games -> all strengths equal at the prior
    assert len(set(out["strengths"].values())) == 1

    # mixed: a beats b decisively in BOTH directions, c ties with everyone
    def judge_mixed(name_a, text_a, name_b, text_b):
        if {text_a, text_b} == {"A", "B"}:
            return (JudgeScores(5, 5, 5, 5) if text_a == "A"
                    else JudgeScores(1, 1, 5, 1))
        return scores  # ties involving C

    out2 = run_tournament({"a": "A", "b": "B", "c": "C"}, judge_mixed,
                          cfg=CampaignConfig(), rng=random.Random(0))
    assert out2["wins"]["a"] == 2 and out2["wins"]["b"] == 0
    assert out2["wins"]["c"] == 0
    assert out2["discarded_games"] == 4  # a-c and b-c, both directions
    assert out2["ranking"][0] == "a"


# ---------------------------------------------------------------------------
# D-9: paper constants are pinned (REPRODUCTION.md fidelity claims)
# ---------------------------------------------------------------------------

def test_paper_constants_are_pinned():
    cfg = CampaignConfig()
    assert cfg.max_concurrent_sessions == 58
    assert cfg.sandbox_cpus == 60
    assert cfg.sandbox_mem_gib == 192
    assert cfg.max_revisions == 10 and cfg.max_gate_failures == 10
    assert cfg.judge_weights == {"impact": 0.35, "novelty": 0.30,
                                 "soundness": 0.25, "actionability": 0.10}
    assert cfg.soundness_auto_lose == 2
    assert cfg.benchmark_attempts_per_cell == 100
    assert cfg.benchmark_max_output_tokens == 1_000_000
    assert cfg.benchmark_levels == ("L1", "L2", "L3", "L4", "L5")


# ---------------------------------------------------------------------------
# C1: the L-ladder is enforced, not data-only
# ---------------------------------------------------------------------------

def test_disallowed_tools_per_level(tmp_path: Path):
    from benchmark.levels import build_environment, disallowed_for

    inputs = tmp_path / "inputs"
    from benchmark.levels import write_synthetic_inputs
    write_synthetic_inputs(inputs, n_loci=4)
    for level, expect in (("L1", True), ("L2", True), ("L3", True),
                          ("L4", True), ("L5", False)):
        env = build_environment(level, inputs)
        dis = disallowed_for(env)
        assert bool(dis) == expect, (level, dis)
        if level == "L3":
            assert set(dis) == {"WebFetch", "WebSearch"}
    l1 = build_environment("L1", inputs)
    assert "Bash" in disallowed_for(l1)


def test_materialize_copies_l3_files_into_the_run_dir(tmp_path: Path):
    # S1 C1: L3+ attempts ran against an EMPTY run dir while the prompt
    # advertised loci/, gene_calls.tsv, pfam_matches.tsv
    from benchmark.levels import build_environment, materialize, write_synthetic_inputs

    inputs = write_synthetic_inputs(tmp_path / "inputs", n_loci=3)
    env = build_environment("L3", inputs)
    run_dir = materialize(env, tmp_path / "run001")
    assert (run_dir / "gene_calls.tsv").exists()
    assert (run_dir / "loci" / "locus_001.fasta").exists()
    assert not (run_dir / "structures").exists()  # L3 has no structures
    env4 = build_environment("L4", inputs)
    run4 = materialize(env4, tmp_path / "run002")
    assert (run4 / "structures" / "locus_001_orf1.pdb").exists()


def test_build_command_carries_the_tool_gate():
    from artharness.runner.base import SessionSpec, build_command

    spec = SessionSpec(role="worker", task_id="b/x/L1/1", system_prompt="s",
                       user_prompt="u", workdir=Path("."),
                       disallowed_tools=["Bash", "WebFetch"])
    cmd = build_command(spec, "model-x")
    assert "--disallowedTools" in cmd
    assert cmd[cmd.index("--disallowedTools") + 1] == "Bash,WebFetch"
    plain = build_command(SessionSpec(role="worker", task_id="t", system_prompt="s",
                                      user_prompt="u", workdir=Path(".")), "m")
    assert "--disallowedTools" not in plain


# ---------------------------------------------------------------------------
# C2: hedged findings must not credit a rubric claim
# ---------------------------------------------------------------------------

def test_hedged_finding_does_not_credit():
    from benchmark.rubric import Finding, Submission, default_rubric, grade_submission

    rubric = default_rubric()
    hedged = Submission(findings=[
        Finding(claim="the 5-prime non-coding tract contains a tandem repeat array",
                evidence="one of several possibilities remains; could be an artifact"),
    ])
    assert grade_submission(hedged, rubric).recognized_repeat_array is False
    assert grade_submission(hedged, rubric).score == 0
    firm = Submission(findings=[
        Finding(claim="the 5-prime non-coding tract contains a tandem repeat array",
                evidence="14-mer array, 12 copies, shuffle-controlled"),
    ])
    assert grade_submission(firm, rubric).recognized_repeat_array is True


def test_rubric_id_in_free_text_is_not_screened_out():
    from benchmark.rubric import (
        Finding,
        Submission,
        default_rubric,
        grade_submission,
        screen_findings,
    )

    rubric = default_rubric()
    f = Finding(claim="nct-1 claim is supported",
                evidence="the 5-prime non-coding tract contains a tandem repeat "
                         "array upstream of the RT; see Fig 3 alignment")
    assert f in screen_findings(Submission(findings=[f]), rubric)
    out = grade_submission(Submission(findings=[f]), rubric)
    assert out.asserted["nct-1"] is True


# ---------------------------------------------------------------------------
# A6/D-8: curated.md lands in the record; supervisor prompt carries KB context
# ---------------------------------------------------------------------------

def test_curator_writes_curated_md(tmp_path: Path):
    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import STAGES, Orchestrator
    from artharness.records import RecordStore
    from artharness.runner.base import ScriptedBackend

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    kb = KnowledgeBase(root / "kb")
    ledger = SessionLedger(root / "ledger.jsonl")
    orch = Orchestrator(CampaignConfig(), store, kb, ledger, ScriptedBackend(),
                        {s: (lambda s: True) for s in STAGES})
    orch.run_stage_chain({"1_input_assembly": ["task"]})
    orch.run()
    rec = store.list_tasks()[0]
    assert (store.records / rec.task_id / "curated.md").exists()


def test_supervisor_prompt_contains_kb_context(tmp_path: Path):
    # S1 D: KB entries reached worker prompts only; the paper's supervisors
    # read the shared knowledge too. Deleting the worker injection must FAIL
    # this test (prompt-content assertions, not backend behavior).
    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import Orchestrator
    from artharness.records import RecordStore
    from artharness.runner.base import ScriptedBackend

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    kb = KnowledgeBase(root / "kb")
    ledger = SessionLedger(root / "ledger.jsonl")
    backend = ScriptedBackend()
    orch = Orchestrator(CampaignConfig(), store, kb, ledger, backend,
                        {s: (lambda s: True) for s in
                         ("1_input_assembly", "2_database_sweep",
                          "3_rt_classification", "4_neighborhood_census",
                          "5_deep_dives")})
    brief = "Survey the Madawaska upstream region for tandem repeats"
    orch.run_stage_chain({"1_input_assembly": [brief]})
    rec = store.list_tasks()[0]
    # the entry exists before BOTH roles run (a later entry cannot be retrieved)
    kb.add("t0000", "prior finding", "Madawaska upstream region: array-like run noted")
    orch.roles.worker(rec)
    orch.roles.supervisor(rec)
    supervisor_prompts = [c.user_prompt for c in backend.calls
                          if c.role == "supervisor"]
    assert any("Madawaska" in p for p in supervisor_prompts)
    worker_prompts = [c.user_prompt for c in backend.calls if c.role == "worker"]
    assert any("Madawaska" in p for p in worker_prompts)


def test_campaign_record_carries_brief_provenance(tmp_path: Path):
    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import STAGES, Orchestrator
    from artharness.records import RecordStore

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    orch = Orchestrator(CampaignConfig(), store, KnowledgeBase(root / "kb"),
                        SessionLedger(root / "ledger.jsonl"), ScriptedBackend(),
                        {s: (lambda s: True) for s in STAGES})
    orch.run_stage_chain({"1_input_assembly": ["task"]})
    text = (root / "campaign.md").read_text()
    assert "reconstruction-v0 (6 anchors), not verbatim" in text
    assert "not the" in text and "verbatim Anthropic brief" in text


# ---------------------------------------------------------------------------
# S3 finding 3: a failed revision check must re-dispatch the WORKER, not hand
# an empty task to an accepting supervisor (curator then crashed on FileNotFound)
# ---------------------------------------------------------------------------

class ReviseOnceThenAccept(ScriptedBackend):
    def __init__(self):
        super().__init__()
        self.worker_calls = 0
        self.sup_calls = 0

    def run(self, spec):
        out = ScriptedBackend.run(self, spec)
        if spec.role == "worker":
            self.worker_calls += 1
            if self.worker_calls == 2:  # one empty revision
                (spec.workdir / "summary.md").unlink(missing_ok=True)
        if spec.role == "supervisor":
            self.sup_calls += 1
            if self.sup_calls == 1:
                out.verdict = "revise"  # forces the empty revision pass
        return out


def test_empty_revision_recovers_via_redispatch_not_curator_crash(tmp_path: Path):
    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import STAGES, Orchestrator
    from artharness.records import RecordStore

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    orch = Orchestrator(CampaignConfig(), store, KnowledgeBase(root / "kb"),
                        SessionLedger(root / "ledger.jsonl"),
                        ReviseOnceThenAccept(),
                        {s: (lambda s: True) for s in STAGES})
    orch.run_stage_chain({"1_input_assembly": ["task"]})
    orch.run()
    rec = store.get("t0001")
    # pre-fix: supervisor reviewed the summary-less task, accepted, and the
    # curator raised FileNotFoundError mid-campaign
    assert rec.status is TaskStatus.CURATED
    assert rec.gate_failures == 1
    assert orch.report.completed == 1


# ---------------------------------------------------------------------------
# S3 finding 1: tasks_total must not double-count (queue is a subset of store)
# ---------------------------------------------------------------------------

def test_tasks_total_counts_each_task_once(tmp_path: Path):
    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import STAGES, Orchestrator
    from artharness.records import RecordStore

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    orch = Orchestrator(CampaignConfig(), store, KnowledgeBase(root / "kb"),
                        SessionLedger(root / "ledger.jsonl"), ScriptedBackend(),
                        {s: (lambda s: True) for s in STAGES})
    orch.run_stage_chain({STAGES[0]: ["a", "b"], STAGES[1]: ["c"]})
    report = orch.run()
    assert report.tasks_total == 3
    assert report.completed == 3


# ---------------------------------------------------------------------------
# S3 finding 2: the report gate is an explicit state check, not a StrEnum `<`
# ---------------------------------------------------------------------------

def test_file_report_refuses_pre_curation_states(tmp_path: Path):
    import pytest

    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import STAGES, Orchestrator
    from artharness.records import RecordStore

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    orch = Orchestrator(CampaignConfig(), store, KnowledgeBase(root / "kb"),
                        SessionLedger(root / "ledger.jsonl"), ScriptedBackend(),
                        {s: (lambda s: True) for s in STAGES})
    orch._new_task("1_input_assembly", "a", TaskOrigin.SEED)
    # OPEN (never dispatched): the old alphabetical `<` happened to refuse this;
    # EXECUTED (worker done, no curator) MUST also be refused and was not
    # reachable through the old comparison
    with pytest.raises(ValueError, match="curator"):
        orch.file_report("t0001", "# early report")
    assert not (store.records / "t0001" / "report.md").exists()


# ---------------------------------------------------------------------------
# S3 finding 6: a skipped copy may come FIRST in the chain
# ---------------------------------------------------------------------------

def test_chains_allow_skip_as_first_gap():
    assert any(len(c) == 3 for c in _chains([0, 400, 600]))
    full = _chains([0, 400, 600, 800, 1000])
    assert any(len(c) == 5 for c in full), full  # one skip, not two


# ---------------------------------------------------------------------------
# S4 finding 1: follow-ups from a FAILED worker pass must not enter triage
# ---------------------------------------------------------------------------

class FailingWorkerWithFollowups(ScriptedBackend):
    """Every worker pass proposes a follow-up and produces no summary."""

    def run(self, spec):
        out = ScriptedBackend.run(self, spec)
        if spec.role == "worker":
            out.proposed_followups = ["new task please"]
            (spec.workdir / "summary.md").unlink(missing_ok=True)
        return out


def test_failed_passes_do_not_spawn_followup_tasks(tmp_path: Path):
    from artharness.accounting import SessionLedger
    from artharness.knowledge import KnowledgeBase
    from artharness.orchestrator import STAGES, Orchestrator
    from artharness.records import RecordStore

    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    orch = Orchestrator(CampaignConfig(max_gate_failures=2, max_tasks_total=30),
                        store, KnowledgeBase(root / "kb"),
                        SessionLedger(root / "ledger.jsonl"),
                        FailingWorkerWithFollowups(),
                        {s: (lambda s: True) for s in STAGES})
    orch.run_stage_chain({"1_input_assembly": ["seed task"]})
    orch.run()
    # pre-fix: every retried pass re-proposed, so the budget filled with
    # follow-up tasks (30 tasks / 29 follow-ups on grok's probe)
    assert len(store.list_tasks()) == 1
    assert orch.report.follow_ups == 0
    assert store.get("t0001").status is TaskStatus.STALLED


# ---------------------------------------------------------------------------
# S4 finding 2: chains must not jump over DETECTED copies
# ---------------------------------------------------------------------------

def _index_jumps_ok(chain: list[int], positions: list[int]) -> bool:
    """At most one index skip per chain, and each skip is a single position."""
    idx = [positions.index(p) for p in chain]
    skips = sum(b - a - 1 for a, b in zip(idx, idx[1:], strict=False))
    return skips <= 1 and all(0 <= b - a <= 2 for a, b in zip(idx, idx[1:],
                                                             strict=False))


def test_chains_never_step_over_detected_copies():
    for positions in ([0, 180, 200, 400, 600, 800],
                      [0, 200, 250, 400, 600, 800]):
        for chain in _chains(positions):
            assert _index_jumps_ok(chain, positions), (chain, positions)
    # grok's probes: the longest legal chains exclude the un-fittable prefix
    assert max(len(c) for c in _chains([0, 180, 200, 400, 600, 800])) == 4
    assert max(len(c) for c in _chains([0, 200, 250, 400, 600, 800])) == 3
    # the S3 no.6 behaviors are preserved
    assert any(len(c) == 3 for c in _chains([0, 400, 600]))
    assert any(len(c) == 5 for c in _chains([0, 400, 600, 800, 1000]))
    assert not any(len(c) == 5 for c in _chains([0, 200, 400, 800, 1200]))
