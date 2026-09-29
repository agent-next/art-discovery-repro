"""End-to-end orchestrator behavior against the offline scripted backend."""

from pathlib import Path

import pytest

from artharness.accounting import SessionLedger, SessionResult
from artharness.config import CampaignConfig
from artharness.knowledge import KnowledgeBase
from artharness.orchestrator import STAGES, Orchestrator
from artharness.records import RecordStore, TaskOrigin, TaskStatus
from artharness.runner.base import BackendOutput, ScriptedBackend


def _sr(role, task_id):
    return SessionResult(role=role, task_id=task_id, duration_s=1.0,
                         input_tokens_uncached=10, output_tokens=10,
                         cache_write_tokens=10)


class AlwaysRevise(ScriptedBackend):
    """Workers submit fine; supervisors never accept (tests the stall path)."""

    def run(self, spec):
        out = super().run(spec)
        if spec.role == "supervisor":
            out.verdict = "revise"
            out.verdict_notes = "again"
        return out


def make_orch(tmp_path: Path, backend=None, gates=None, cfg=None, triage=None):
    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=False)
    kb = KnowledgeBase(root / "kb")
    ledger = SessionLedger(root / "ledger.jsonl")
    backend = backend or ScriptedBackend()
    gates = gates or {s: (lambda s: True) for s in STAGES}
    orch = Orchestrator(cfg or CampaignConfig(), store, kb, ledger, backend, gates,
                        triage=triage)
    return orch, store, ledger, backend


def test_stage_chain_gates_block_later_stages(tmp_path: Path):
    unlocked = {s: False for s in STAGES}

    def gate_for(stage):
        return lambda s: unlocked[s]

    gates = {s: gate_for(s) for s in STAGES}
    orch, store, _, _ = make_orch(tmp_path, gates=gates)
    unlocked[STAGES[0]] = True
    with pytest.raises(RuntimeError, match="gate for"):
        orch.run_stage_chain({STAGES[0]: ["assemble"], STAGES[2]: ["classify"]})
    ids = [r.task_id for r in store.list_tasks()]
    assert ids == ["t0001"]  # stage-3 task never opened while stage-2 gate is shut


def test_stage_gates_evaluate_after_stage_work(tmp_path: Path):
    """Power test (paper p.28): a stage's gate sees that stage's finished tasks, and
    no later-stage task exists until it passes. The previous implementation opened
    every stage before dispatching anything, so each gate ran against zero work and
    the final stage's gate never ran."""
    events: list[tuple[str, dict[str, str | None]]] = []
    holder: dict = {}

    def gate_for(stage):
        def gate(_):
            snap = {r.task_id: (r.stage, r.status.value)
                    for r in holder["store"].list_tasks()}
            events.append((stage, snap))
            return True
        return gate

    gates = {s: gate_for(s) for s in STAGES}
    orch, store, _, _ = make_orch(tmp_path, gates=gates)
    holder["store"] = store
    orch.run_stage_chain({STAGES[0]: ["assemble"], STAGES[1]: ["sweep"],
                          STAGES[4]: ["dive"]})

    assert [e[0] for e in events] == list(STAGES)
    for stage, snap in events:
        idx = STAGES.index(stage)
        done = [st for (stg, st) in snap.values() if STAGES.index(stg) <= idx]
        assert done and all(st == TaskStatus.CURATED.value for st in done)
        assert all(STAGES.index(stg) <= idx for (stg, _) in snap.values())


def test_full_happy_path_with_curation(tmp_path: Path):
    orch, store, ledger, _ = make_orch(tmp_path)
    orch.run_stage_chain({STAGES[0]: ["Assemble input data"],
                          STAGES[4]: ["Deep dive A"]})
    report = orch.run()
    assert report.completed == 2 and report.stalled == 0
    for rec in store.list_tasks():
        assert rec.status in (TaskStatus.CURATED, TaskStatus.COMPLETED)
        assert (store.records / rec.task_id / "summary.md").exists()
        assert list((store.root / "kb" / "entries").glob("*.md"))
    roles = {s.role for s in ledger.sessions()}
    assert {"worker", "supervisor", "curator"} <= roles


def test_revision_then_accept(tmp_path: Path):
    class ReviseOnce(ScriptedBackend):
        """First supervisor call revises, later ones accept (the paper's mode:
        49/119 tasks revised at least once, then accepted)."""

        def __init__(self):
            super().__init__()
            self.revised_once = False

        def run(self, spec):
            out = super().run(spec)
            if spec.role == "supervisor" and not self.revised_once:
                self.revised_once = True
                out.verdict = "revise"
                out.verdict_notes = "recount repeats"
            return out

    orch, store, _, _ = make_orch(tmp_path, backend=ReviseOnce())
    orch.run_stage_chain({STAGES[0]: ["single task"]})
    report = orch.run()
    assert report.revised == 1 and report.completed == 1
    rec = store.get("t0001")
    assert rec.revisions == 1
    assert rec.status in (TaskStatus.CURATED, TaskStatus.COMPLETED)
    assert "revise" in store.read_text("t0001", "verdict-r1.md")
    assert "accept" in store.read_text("t0001", "verdict.md")


def test_stall_after_max_revisions(tmp_path: Path):
    orch, store, _, _ = make_orch(tmp_path, backend=AlwaysRevise(),
                                  cfg=CampaignConfig(max_revisions=3))
    orch.run_stage_chain({STAGES[0]: ["never good enough"]})
    report = orch.run()
    assert report.stalled == 1 and report.completed == 0
    assert store.get("t0001").status is TaskStatus.STALLED


def test_stall_after_failed_completion_checks(tmp_path: Path):
    orch, store, _, _ = make_orch(
        tmp_path, backend=ScriptedBackend(write_outputs=False),
        cfg=CampaignConfig(max_gate_failures=2))
    orch.run_stage_chain({STAGES[0]: ["worker never produces a summary"]})
    report = orch.run()
    assert report.stalled == 1
    assert store.get("t0001").status is TaskStatus.STALLED


def test_followup_triage_accept_and_reject(tmp_path: Path):
    def triage(brief, parent):
        if "genome" in brief:
            return True, ""
        return False, "outside mission scope"

    scenario = {("worker", "t0001"): BackendOutput(
        result=_sr("worker", "t0001"),
        proposed_followups=["study genome neighborhood X", "coffee break task"])}
    orch, store, _, _ = make_orch(tmp_path, backend=ScriptedBackend(scenario),
                                  triage=triage)
    orch.run_stage_chain({STAGES[0]: ["seed"]})
    report = orch.run()
    assert report.rejected_at_triage == 1
    rejections = list((store.records / "t0001").glob("triage-rejection-*.md"))
    assert len(rejections) == 1
    assert "outside mission scope" in rejections[0].read_text()
    followups = [r for r in store.list_tasks() if r.origin is TaskOrigin.FOLLOW_UP]
    assert len(followups) == 1 and followups[0].parent == "t0001"
    assert report.follow_ups == 1


def test_report_editor_gate(tmp_path: Path):
    scenario = {("editor", "t0002"): BackendOutput(
        result=_sr("editor", "t0002"), verdict="no",
        verdict_notes="soundness of claim 3 unestablished")}
    orch, store, _, _ = make_orch(tmp_path, backend=ScriptedBackend(scenario))
    orch.run_stage_chain({STAGES[0]: ["a"], STAGES[1]: ["b"]})
    orch.run()
    filed = orch.file_report("t0001", "# Report A\nsolid")
    rejected = orch.file_report("t0002", "# Report B\nshaky")
    assert filed is not None and filed.exists() and filed.name == "report.md"
    assert rejected is None
    assert "FILE: no" in store.read_text("t0002", "report-review.md")
    assert orch.report.reports_filed == 1
    assert store.get("t0002").status is TaskStatus.CURATED  # unchanged by rejection


class RevisedTaskRepeatsFollowup(ScriptedBackend):
    """Task t0001's worker proposes the same follow-up on every pass and the
    supervisor revises once before accepting."""

    def __init__(self):
        super().__init__()
        self.supervisor_calls = 0

    def run(self, spec):
        out = super().run(spec)
        if spec.role == "worker" and spec.task_id == "t0001":
            out.proposed_followups = ["follow up on X"]
        if spec.role == "supervisor" and spec.task_id == "t0001":
            self.supervisor_calls += 1
            if self.supervisor_calls == 1:
                out.verdict = "revise"
        return out


def test_revision_passes_do_not_duplicate_followups(tmp_path: Path):
    orch, store, _, _ = make_orch(tmp_path, backend=RevisedTaskRepeatsFollowup())
    orch.run_stage_chain({STAGES[0]: ["seed"]})
    tasks = store.list_tasks()
    assert store.get("t0001").revisions == 1  # two worker passes happened
    assert [t.origin for t in tasks].count(TaskOrigin.FOLLOW_UP) == 1
    assert orch.report.follow_ups == 1


def test_launch_session_recorded_and_task_mix(tmp_path: Path):
    # Paper p.30: sessions = the LAUNCH session plus the role sessions
    # (949 = 1 + 414 + 375 + 107 + 52); task mix = 1 assembly + 96 analyses
    # + 22 report writing. The orchestrator used to record neither the
    # launch session nor any task-mix axis. FAILING-FIRST.
    orch, store, ledger, _ = make_orch(tmp_path)
    orch.run_stage_chain({STAGES[0]: ["assemble"], STAGES[1]: ["sweep"],
                          STAGES[2]: ["classify"]})
    roles = ledger.summary()["sessions_by_role"]
    assert roles["launch"] == 1  # exactly one launch session per campaign
    mix = store.task_mix()
    assert mix == {"assembly": 1, "analysis": 2}
    rec = store.list_tasks()[0]
    assert rec.kind.value == "assembly"


def test_revision_cascade_is_bounded_by_followup_dedup(tmp_path: Path):
    # Internal-audit finding 7 (2026-09-28): a revised task re-presenting the
    # SAME follow-up on every pass used to multiply tasks (offline repro hit
    # max_tasks_total=500 with 499 identical briefs). The (parent, brief)
    # dedup in propose_followup must bound the cascade: one seed task that
    # always proposes one identical follow-up and always gets revised must
    # produce exactly ONE follow-up task, never a chain.
    # drive the orchestrator API directly: the cascade is a property of
    # propose_followup + revisions, not of the backend
    orch, store, ledger, _ = make_orch(tmp_path)
    orch.run_stage_chain({STAGES[0]: ["seed brief"]})
    seed = store.list_tasks()[0]
    for _ in range(20):  # 20 revision passes, same follow-up each time
        orch.propose_followup("identical follow-up brief", seed)
    follow_ups = [r for r in store.list_tasks()
                  if r.origin.value == "follow_up"]
    assert len(follow_ups) == 1, [r.brief for r in follow_ups]
    # distinct briefs still open their own tasks
    orch.propose_followup("a different follow-up brief", seed)
    follow_ups = [r for r in store.list_tasks()
                  if r.origin.value == "follow_up"]
    assert len(follow_ups) == 2
