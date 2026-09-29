"""The campaign orchestrator: stage chain, task queue, triage, dispatch loop.

Faithful mechanics from the paper (Methods p.28, "Autonomous research harness"):

- launch stage chain: "A launch agent encoded the five stages of the research brief as
  a chain, with each stage closed by scripted completion checks... No task of a later
  stage was opened until the preceding stage was completed."
- worker/supervisor loop: 49 of 119 tasks revised at least once; stalls recorded
  "after ten revisions" / "after ten failed completion checks".
- triage: "Follow-up tasks proposed by agents during the work entered a triage queue,
  from which the harness released or rejected with a written reason" (10 rejected).
- curation after each completed task; editor review before filing reports.
- termination: "A campaign ended when no task remained that could be dispatched."

Gates are *scripted* callables supplied by the caller (the paper's completion checks
were scripts, not model judgments). The orchestrator never invents scientific
acceptance criteria.
"""

from __future__ import annotations

import threading
import time
from collections import deque
from collections.abc import Callable
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path

from .accounting import SessionLedger, SessionResult
from .config import CampaignConfig
from .knowledge import KnowledgeBase
from .records import RecordStore, TaskOrigin, TaskRecord, TaskStatus
from .roles import Roles

STAGES = (
    "1_input_assembly",
    "2_database_sweep",
    "3_rt_classification",
    "4_neighborhood_census",
    "5_deep_dives",
)

# Gate signature: stage_name -> True when the stage's scripted completion check passes.
Gate = Callable[[str], bool]

# Triage signature: (brief, parent_record) -> (release: bool, reason: str).
# A rejection MUST carry a written reason (paper: "released or rejected with a
# written reason"). Default policy: release everything; budget caps are enforced
# separately in Orchestrator._new_task.
TriagePolicy = Callable[[str, TaskRecord | None], tuple[bool, str]]


@dataclass
class CampaignReport:
    tasks_total: int = 0
    completed: int = 0
    rejected_at_triage: int = 0
    stalled: int = 0
    revised: int = 0
    follow_ups: int = 0
    reports_filed: int = 0

    def render(self) -> str:
        lines = [f"tasks_total={self.tasks_total} completed={self.completed} "
                 f"rejected_at_triage={self.rejected_at_triage} stalled={self.stalled}",
                 f"revised>=1: {self.revised}  follow_ups_opened: {self.follow_ups}  "
                 f"reports_filed: {self.reports_filed}"]
        return "\n".join(lines)


class Orchestrator:
    def __init__(self, cfg: CampaignConfig, store: RecordStore, kb: KnowledgeBase,
                 ledger: SessionLedger, backend, gates: dict[str, Gate],
                 triage: TriagePolicy | None = None, skills_dir: Path | None = None):
        self.cfg = cfg
        self.store = store
        self.kb = kb
        self.ledger = ledger
        self.roles = Roles(backend, store, kb, skills_dir)
        self.gates = gates
        self.triage = triage or (lambda brief, parent: (True, ""))
        self.queue: deque[str] = deque()
        # Guards every write to shared campaign state (task creation, queue, report
        # counters, ledger, triage bookkeeping). Backend sessions run OUTSIDE it.
        self._lock = threading.RLock()
        self._proposed: set[tuple[str, str]] = set()
        self.report = CampaignReport()

    # -- stage chain ---------------------------------------------------------
    def run_stage_chain(self, briefs: dict[str, list[str]],
                        deep_dive_labels: dict[str, list[str]] | None = None) -> None:
        """Execute the research brief as a chain of stages (paper: 5 stage tasks, 16
        deep dives). A stage's tasks are opened only after the previous stage's
        tasks (and every follow-up they spawned) have been dispatched and its
        scripted gate passes; the last stage's gate closes the chain."""
        # brief provenance is part of the versioned record (S1 finding D: the
        # mandate existed only as prose in the brief; nothing enforced it)
        # The launch session (paper p.30: 949 sessions = the launch session
        # plus 414 worker / 375 supervisor / 107 curator / 52 editor) is part
        # of the campaign accounting; record it around the whole chain.
        # NOT-IN-PAPER: token split for the launch session is not observable
        # at this layer — zeros, never fabricated.
        launched_at = time.monotonic()
        try:
            self._run_stage_chain_inner(briefs, deep_dive_labels)
        finally:
            with self._lock:
                self.ledger.record(SessionResult(
                    role="launch", task_id=None,
                    duration_s=time.monotonic() - launched_at,
                    input_tokens_uncached=0, output_tokens=0,
                    cache_write_tokens=0))

    def _run_stage_chain_inner(self, briefs: dict[str, list[str]],
                               deep_dive_labels: dict[str, list[str]] | None) -> None:
        (self.store.root / "campaign.md").write_text(
            f"# campaign record\n\n{self.cfg.brief_provenance}\n\n"
            "This campaign runs a RECONSTRUCTED research brief; it is not the "
            "verbatim Anthropic brief (Supplementary Note 1, unpublished).\n")
        deep_dive_labels = deep_dive_labels or {}
        for i, stage in enumerate(STAGES):
            if i:
                self._require_gate(STAGES[i - 1])
            for brief in briefs.get(stage, []):
                self._new_task(stage, brief, TaskOrigin.SEED)
            for label in deep_dive_labels.get(stage, []):
                self._new_task(stage, f"Deep dive: {label}", TaskOrigin.DEEP_DIVE)
            self._drain()
        self._require_gate(STAGES[-1])
        self.report.tasks_total = len(self.store.list_tasks())

    def _require_gate(self, stage: str) -> None:
        gate = self.gates.get(stage)
        if gate is None:
            raise ValueError(f"stage {stage} has no scripted completion gate")
        if not gate(stage):
            raise RuntimeError(f"gate for {stage} failed; aborting chain")

    # -- task creation / triage ------------------------------------------------
    def _new_task(self, stage: str, brief: str, origin: TaskOrigin,
                  parent: str | None = None, label: str = "") -> TaskRecord | None:
        with self._lock:
            n = len(self.store.list_tasks())
            if n >= self.cfg.max_tasks_total:
                self._reject_uncreated(brief, "task budget exhausted (max_tasks_total)")
                return None
            rec = self.store.create(self.store.next_task_id(n + 1),
                                    brief, stage, origin, parent=parent, label=label)
            self.queue.append(rec.task_id)
            if origin is TaskOrigin.FOLLOW_UP:
                self.report.follow_ups += 1
            return rec

    def _reject_uncreated(self, brief: str, reason: str) -> None:
        # Budget rejections are still "rejected with a written reason"; the reason is
        # recorded in the campaign log (no task dir exists to hold it).
        self.report.rejected_at_triage += 1
        with (self.store.root / "triage-rejections.log").open("a") as log:
            log.write(f"REJECTED: {reason}\nbrief: {brief[:400]}\n\n")

    def propose_followup(self, brief: str, parent: TaskRecord,
                         proposed_by: str = "worker") -> None:
        with self._lock:
            # A revised task re-presents the same follow-ups on every pass; triage each
            # distinct (parent, brief) once so revisions cannot multiply tasks.
            key = (parent.task_id, brief)
            if key in self._proposed:
                return
            self._proposed.add(key)
            ok, reason = self.triage(brief, parent)
            if ok:
                self._new_task(parent.stage, brief, TaskOrigin.FOLLOW_UP, parent=parent.task_id,
                               label=f"followup by {proposed_by}")
            else:
                self.report.rejected_at_triage += 1
                seq = sum(1 for _ in
                          self.store.records.glob(parent.task_id + "/triage-rejection-*.md")) + 1
                self.store.write_text(
                    parent.task_id, f"triage-rejection-{seq}.md",
                    f"REJECTED AT TRIAGE\nreason: {reason}\n\nbrief:\n{brief}\n")
                self.store.commit(f"triage({parent.task_id}): reject follow-up — {reason}")

    # -- dispatch loop ------------------------------------------------------------
    def run(self) -> CampaignReport:
        """Run until the queue is exhausted (paper's termination condition)."""
        self._drain()
        # S3 finding 1: the queue is a SUBSET of the store (run_stage_chain
        # appends every created task to both) — adding them double-counted.
        self.report.tasks_total = len(self.store.list_tasks())
        return self.report

    def _drain(self) -> None:
        workers = max(1, min(self.cfg.dispatch_concurrency,
                             self.cfg.max_concurrent_sessions))
        if workers == 1:
            while self.queue:
                self._dispatch(self.store.get(self.queue.popleft()))
            return
        # Tasks opened while others run (follow-ups) join the pool; a stage is drained
        # only when nothing is queued AND nothing is in flight.
        with ThreadPoolExecutor(max_workers=workers) as pool:
            pending: set = set()
            while True:
                with self._lock:
                    while self.queue:
                        rec = self.store.get(self.queue.popleft())
                        pending.add(pool.submit(self._dispatch, rec))
                if not pending:
                    return
                done, pending = wait(pending, timeout=0.05, return_when=FIRST_COMPLETED)
                for fut in done:
                    fut.result()

    def _dispatch(self, rec: TaskRecord) -> None:
        # One loop per WORKER pass. Every pass — initial or revision — must pass
        # the completion check before the supervisor reviews it (S3 findings 3/4:
        # the old revision path re-entered the supervisor loop with no summary,
        # so an accepting supervisor sent an empty task to the curator and the
        # completion-check stall mode was unreachable after any revision).
        while True:
            out = self.roles.worker(rec)
            self._record(out.result)

            summary = self.store.records / rec.task_id / "summary.md"
            if not summary.exists():
                # A worker pass with no summary counts as a failed completion check
                # (paper stall mode 2: "after ten failed completion checks").
                rec.gate_failures += 1
                rec.status = TaskStatus.OPEN
                self.store.update(rec, f"task({rec.task_id}): completion check failed "
                                       f"({rec.gate_failures}/{self.cfg.max_gate_failures})")
                if rec.gate_failures >= self.cfg.max_gate_failures:
                    self._stall(rec, "completion checks")
                    return
                # S4 finding 1: follow-ups from a FAILED pass must not enter
                # triage — with retry loops a failing worker would re-propose
                # every pass and explode the task budget. Triage happens only
                # for a pass that passed the check (below).
                continue  # re-dispatch the worker

            # worker-proposed follow-ups enter triage (valid passes only)
            for fu in out.proposed_followups:
                self.propose_followup(fu, rec, proposed_by="worker")

            rec.status = TaskStatus.EXECUTED
            self.store.update(rec, f"task({rec.task_id}): worker submitted summary")

            # supervisor review. Only an explicit "accept" accepts; a missing or
            # unparseable verdict fails toward revision (grok review 2026-09-24:
            # the old code treated None as accept, so live backends that omitted
            # the verdict could never trigger revisions or stalls).
            revise = True
            while revise:
                sout = self.roles.supervisor(rec)
                self._record(sout.result)
                # supervisor-authored follow-up briefs enter triage too (paper
                # A1: the t0010 supervisor wrote the t0062 brief that led to ART)
                for fu in sout.proposed_followups:
                    self.propose_followup(fu, rec, proposed_by="supervisor")
                if sout.verdict != "accept":
                    rec.revisions += 1
                    rec.status = TaskStatus.REVISING
                    self.store.write_text(
                        rec.task_id, f"verdict-r{rec.revisions}.md",
                        f"VERDICT: revise ({sout.verdict or 'no verdict parsed'})"
                        f"\n\n{sout.verdict_notes or ''}\n")
                    self.store.update(rec, f"task({rec.task_id}): supervisor revise "
                                           f"({rec.revisions}/{self.cfg.max_revisions})")
                    if rec.revisions >= self.cfg.max_revisions:
                        self._stall(rec, "revisions")
                        return
                    # The pass-1 summary is a stale output once "revise" is
                    # issued: withdraw it so the revision pass must pass the
                    # same completion check instead of silently re-presenting
                    # the old summary.
                    summary.unlink(missing_ok=True)
                    break  # back to the worker pass (outer loop)
                revise = False

            if not revise:
                self.store.write_text(rec.task_id, "verdict.md",
                                      f"VERDICT: accept\n\n{sout.verdict_notes or ''}\n")
                rec.status = TaskStatus.ACCEPTED
                self.store.update(rec, f"task({rec.task_id}): supervisor accepted")
                if rec.revisions >= 1:
                    self._bump("revised")

                # curator enters findings into the shared knowledge base
                cout = self.roles.curator(rec)
                self._record(cout.result)
                rec.status = TaskStatus.CURATED
                self.store.update(rec, f"task({rec.task_id}): curated")

                self._bump("completed")
                return

    def _bump(self, counter: str) -> None:
        with self._lock:
            setattr(self.report, counter, getattr(self.report, counter) + 1)

    def _record(self, result) -> None:
        with self._lock:
            self.ledger.record(result)

    def _stall(self, rec: TaskRecord, mode: str) -> None:
        rec.status = TaskStatus.STALLED
        self.store.update(rec, f"task({rec.task_id}): STALLED after {mode}")
        if rec.revisions >= 1:
            # paper counts a task "revised >= 1x" whenever revisions happened,
            # including tasks that later stalled (49/119 figure)
            self._bump("revised")
        self._bump("stalled")

    # -- reports -----------------------------------------------------------------
    def file_report(self, task_id: str, report_text: str) -> Path | None:
        """Write a draft report and route it through editor review (52 editor sessions
        for 19 reports in the paper). Editor-rejected reports are NOT filed."""
        rec = self.store.get(task_id)
        # S3 finding 2: TaskStatus is a StrEnum — `<` compares the STRING values
        # alphabetically, not the lifecycle order, so the old guard accepted
        # ACCEPTED (never curated) and refused EXECUTED... by accident. Compare
        # states explicitly.
        if rec.status is not TaskStatus.CURATED:
            raise ValueError(f"task {task_id} not ready for a report "
                             f"(status {rec.status}); the curator step must "
                             f"complete first")
        draft = self.store.write_text(task_id, "report-draft.md", report_text)
        eout = self.roles.editor(task_id, draft)
        self._record(eout.result)
        if eout.verdict != "yes":  # only an explicit yes files; None fails safe
            self.store.write_text(task_id, "report-review.md",
                                  f"FILE: no\n\n{eout.verdict_notes or ''}\n")
            self.store.commit(f"report({task_id}): editor rejected")
            return None
        final = self.store.write_text(task_id, "report.md", report_text)
        self.store.commit(f"report({task_id}): filed after editor review")
        rec.status = TaskStatus.COMPLETED
        self.store.update(rec, f"task({task_id}): report filed")
        self._bump("reports_filed")
        return final
