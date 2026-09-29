"""Threaded dispatch (paper Methods p.28: up to 58 concurrent sessions).

Each test pins a property a sequential loop cannot satisfy or a racy one would break:
real overlap, the concurrency cap, and integrity of shared campaign state.
"""

import subprocess
import threading
import time
from pathlib import Path

import pytest

from artharness.accounting import SessionLedger
from artharness.config import CampaignConfig
from artharness.knowledge import KnowledgeBase
from artharness.orchestrator import STAGES, Orchestrator
from artharness.records import RecordStore, TaskOrigin, TaskStatus
from artharness.runner.base import ScriptedBackend


def make_orch(tmp_path: Path, backend, cfg: CampaignConfig, use_git: bool = False):
    root = tmp_path / "campaign"
    root.mkdir()
    store = RecordStore(root, use_git=use_git)
    orch = Orchestrator(cfg, store, KnowledgeBase(root / "kb"),
                        SessionLedger(root / "ledger.jsonl"), backend,
                        {s: (lambda s: True) for s in STAGES})
    return orch, store


class BarrierBackend(ScriptedBackend):
    """Every worker waits for `parties` workers to be in flight at once."""

    def __init__(self, parties: int, timeout: float):
        super().__init__()
        self.barrier = threading.Barrier(parties, timeout=timeout)

    def run(self, spec):
        if spec.role == "worker":
            self.barrier.wait()
        return super().run(spec)


class InflightBackend(ScriptedBackend):
    def __init__(self, delay: float = 0.05):
        super().__init__()
        self.delay = delay
        self.now = 0
        self.peak = 0
        self._mu = threading.Lock()

    def run(self, spec):
        if spec.role == "worker":
            with self._mu:
                self.now += 1
                self.peak = max(self.peak, self.now)
            time.sleep(self.delay)
            with self._mu:
                self.now -= 1
        return super().run(spec)


def test_tasks_overlap_when_concurrency_allows(tmp_path: Path):
    orch, store = make_orch(tmp_path, BarrierBackend(4, timeout=5),
                            CampaignConfig(dispatch_concurrency=4))
    orch.run_stage_chain({STAGES[0]: [f"seed {i}" for i in range(4)]})
    assert orch.report.completed == 4  # the barrier only opens if 4 ran together


def test_sequential_default_cannot_overlap(tmp_path: Path):
    """Power partner of the test above: with the default (1) the same barrier breaks."""
    orch, _ = make_orch(tmp_path, BarrierBackend(4, timeout=0.3), CampaignConfig())
    with pytest.raises(threading.BrokenBarrierError):
        orch.run_stage_chain({STAGES[0]: [f"seed {i}" for i in range(4)]})


def test_concurrency_is_capped(tmp_path: Path):
    backend = InflightBackend()
    orch, _ = make_orch(tmp_path, backend, CampaignConfig(dispatch_concurrency=2))
    orch.run_stage_chain({STAGES[0]: [f"seed {i}" for i in range(8)]})
    assert backend.peak == 2  # used the allowance, never exceeded it

    (tmp_path / "b").mkdir()
    backend = InflightBackend()
    orch, _ = make_orch(tmp_path / "b", backend,
                        CampaignConfig(dispatch_concurrency=10, max_concurrent_sessions=3))
    orch.run_stage_chain({STAGES[0]: [f"seed {i}" for i in range(9)]})
    assert backend.peak == 3  # dispatch_concurrency cannot exceed the paper cap


class FollowUpBackend(ScriptedBackend):
    """Each seed task proposes one uniquely-worded follow-up."""

    def run(self, spec):
        out = super().run(spec)
        if spec.role == "worker" and (spec.workdir / "brief.md").read_text().startswith("seed"):
            out.proposed_followups = [f"child of {spec.task_id}"]
            time.sleep(0.01)
        return out


def test_shared_state_stays_consistent_under_concurrency(tmp_path: Path):
    n = 24
    orch, store = make_orch(tmp_path, FollowUpBackend(),
                            CampaignConfig(dispatch_concurrency=8))
    orch.run_stage_chain({STAGES[0]: [f"seed {i}" for i in range(n)]})
    tasks = store.list_tasks()
    ids = [t.task_id for t in tasks]
    assert len(ids) == len(set(ids)) == 2 * n
    assert all(t.status is TaskStatus.CURATED for t in tasks)
    assert [t.origin for t in tasks].count(TaskOrigin.FOLLOW_UP) == n
    sessions = orch.ledger.sessions()
    # worker + supervisor + curator per task, PLUS the one launch session
    # the chain records (paper p.30 accounting; batch2 2026-09-29)
    assert len(sessions) == 2 * n * 3 + 1
    assert sum(1 for s in sessions if s.role == "launch") == 1
    r = orch.report
    assert (r.completed, r.tasks_total, r.follow_ups, r.stalled) == (2 * n, 2 * n, n, 0)


def test_task_budget_is_not_overshot_under_concurrency(tmp_path: Path):
    cfg = CampaignConfig(dispatch_concurrency=8, max_tasks_total=30)
    orch, store = make_orch(tmp_path, FollowUpBackend(), cfg)
    orch.run_stage_chain({STAGES[0]: [f"seed {i}" for i in range(24)]})
    assert len(store.list_tasks()) == 30  # 24 seeds + only 6 of the 24 follow-ups fit
    assert orch.report.rejected_at_triage == 18


def test_git_audit_trail_loses_no_transition(tmp_path: Path, monkeypatch, capsys):
    for k, v in {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
                 "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"}.items():
        monkeypatch.setenv(k, v)
    root = tmp_path / "campaign"
    root.mkdir()
    subprocess.run(["git", "init", "-q", "-b", "art-test"], cwd=root, check=True)
    store = RecordStore(root, use_git=True)
    assert store.use_git
    orch = Orchestrator(CampaignConfig(dispatch_concurrency=8), store,
                        KnowledgeBase(root / "kb"), SessionLedger(root / "ledger.jsonl"),
                        ScriptedBackend(), {s: (lambda s: True) for s in STAGES})
    orch.run_stage_chain({STAGES[0]: [f"seed {i}" for i in range(16)]})
    assert "git commit failed" not in capsys.readouterr().err
    dirty = subprocess.run(["git", "status", "--porcelain", "records"], cwd=root,
                           capture_output=True, text=True, check=True).stdout
    assert dirty.strip() == ""  # every transition landed in a commit
