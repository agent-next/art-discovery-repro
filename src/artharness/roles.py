"""Role drivers: worker / supervisor / curator / editor.

Paper (Results p.1-2, Methods p.28):
- worker: "proposes a plan and executes it", submits "a written summary of the results
  with the code and data files"; may propose follow-up tasks.
- supervisor: "reviewed the plan, the summary, and the edited files, and either
  accepted them or returned them for revision" (49/119 revised at least once);
  supervisors wrote the briefs of follow-up tasks (the ART chain t0010 -> t0062).
- curator: read brief + summary, "entered the findings into a shared knowledge base".
- editor: "Reports were reviewed by an editor agent before filing" (52 editor
  sessions for 19 reports).

Each driver builds the prompts and delegates cognition to the backend; all side
effects land in the RecordStore. Follow-up task briefs proposed by workers are
verbatim queued to triage; briefs written by supervisors land with
origin=FOLLOW_UP and parent set (paper's t0062 pattern).
"""

from __future__ import annotations

from pathlib import Path

from .knowledge import KnowledgeBase
from .records import RecordStore, TaskRecord
from .runner.base import BackendOutput, SessionBackend, SessionSpec  # noqa: F401

WORKER_SYSTEM = """\
You are the worker agent of a genome-mining research campaign.
Write a concise plan, execute it with the tools in your sandbox, then submit a
written summary with the code and data files you produced (under artifacts/).
Method guides live in the skills directory; consult the relevant one before coding.
If observations suggest a new question worth a dedicated task, PROPOSE_FOLLOWUP with
a one-paragraph task brief; do not expand the current task's scope instead.
"""

SUPERVISOR_SYSTEM = """\
You are the supervisor agent of a genome-mining research campaign.
Review the worker's plan, summary, and edited files. Return VERDICT: accept, or
VERDICT: revise with concrete required changes. Check quantitative claims against
the artifacts. If your review suggests a new task, write its brief after the verdict
under PROPOSED_TASK_BRIEF.
"""

CURATOR_SYSTEM = """\
You are the curator agent of a genome-mining research campaign.
Read the task brief and the worker summary; write a knowledge-base entry stating the
findings, the evidence, and the task ids they came from. Be terse and factual.
"""

EDITOR_SYSTEM = """\
You are the editor agent of a genome-mining research campaign.
Review the draft report for soundness, novelty claims, and unsupported statements.
Return FILE: yes or FILE: no with required edits.
"""


class Roles:
    def __init__(self, backend, store: RecordStore, kb: KnowledgeBase, skills_dir: Path | None):
        self.backend = backend
        self.store = store
        self.kb = kb
        self.skills_dir = skills_dir

    # -- worker ------------------------------------------------------------
    def worker(self, rec: TaskRecord) -> BackendOutput:
        brief = self.store.read_text(rec.task_id, "brief.md")
        kb_block = self.kb.context_block(brief)
        # a revision pass must see why the supervisor returned the last one
        verdicts = sorted((self.store.records / rec.task_id).glob("verdict-r*.md"),
                          key=lambda p: int(p.stem.rsplit("-r", 1)[1]))
        feedback = (f"## supervisor feedback to address (revision {rec.revisions})\n"
                    f"{verdicts[-1].read_text()}\n") if verdicts else ""
        prompt = (
            f"# Task {rec.task_id} ({rec.stage})\n\n{brief}\n\n{feedback}{kb_block}"
            "Submit plan.md, summary.md, and artifacts. Propose follow-ups if warranted."
        )
        out = self.backend.run(SessionSpec(
            role="worker", task_id=rec.task_id, system_prompt=WORKER_SYSTEM,
            user_prompt=prompt, workdir=self.store.records / rec.task_id,
            skills_dir=self.skills_dir,
        ))
        return out

    # -- supervisor ----------------------------------------------------------
    def supervisor(self, rec: TaskRecord) -> BackendOutput:
        parts = [f"# Task {rec.task_id} review\n",
                 "## brief\n", self.store.read_text(rec.task_id, "brief.md"), "\n"]
        if (self.store.records / rec.task_id / "plan.md").exists():
            parts += ["## plan\n", self.store.read_text(rec.task_id, "plan.md"), "\n"]
        if (self.store.records / rec.task_id / "summary.md").exists():
            parts += ["## summary\n", self.store.read_text(rec.task_id, "summary.md"), "\n"]
        artifacts = sorted(p.name for p in self.store.artifact_dir(rec.task_id).iterdir())
        parts += ["## artifacts\n", "\n".join(artifacts) or "(none)"]
        # paper: the shared knowledge base is "readable by every agent" —
        # supervisors judged against what the campaign already knew (S1 finding
        # D: only workers saw kb.context_block)
        context = self.kb.context_block(self.store.read_text(rec.task_id, "brief.md"))
        if context:
            parts += ["\n## relevant knowledge base entries\n", context, "\n"]
        out = self.backend.run(SessionSpec(
            role="supervisor", task_id=rec.task_id, system_prompt=SUPERVISOR_SYSTEM,
            user_prompt="".join(parts), workdir=self.store.records / rec.task_id,
        ))
        return out

    # -- curator ---------------------------------------------------------------
    def curator(self, rec: TaskRecord) -> BackendOutput:
        prompt = (
            f"# Task {rec.task_id}\n\n## brief\n{self.store.read_text(rec.task_id, 'brief.md')}"
            f"\n\n## worker summary\n{self.store.read_text(rec.task_id, 'summary.md')}"
        )
        out = self.backend.run(SessionSpec(
            role="curator", task_id=rec.task_id, system_prompt=CURATOR_SYSTEM,
            user_prompt=prompt, workdir=self.store.records / rec.task_id,
        ))
        # live backends return the entry prose in `text`; offline backends fall
        # back to a provenance stub (grok round-2 finding 3)
        entry = out.text or self._synthesized_entry(out, rec)
        title = f"findings from {rec.task_id} ({rec.label or rec.stage})"
        self.kb.add(rec.task_id, title, entry)
        # the entry also lands in the shared task record (records.py layout
        # documents curated.md; the KB alone is not the versioned record)
        self.store.write_text(rec.task_id, "curated.md",
                              f"# {title}\n\n{entry}\n")
        return out

    # -- editor ------------------------------------------------------------------
    def editor(self, task_id: str, report_path: Path) -> BackendOutput:
        out = self.backend.run(SessionSpec(
            role="editor", task_id=task_id, system_prompt=EDITOR_SYSTEM,
            user_prompt=report_path.read_text(), workdir=self.store.root,
        ))
        return out

    @staticmethod
    def _synthesized_entry(out: BackendOutput, rec: TaskRecord) -> str:
        """Offline backends produce no prose; the entry then records the flow itself.

        Live backends should return the entry text in ``verdict_notes``.
        """
        if out.verdict_notes:
            return out.verdict_notes
        return (
            f"Task {rec.task_id} ({rec.stage}, {rec.origin.value}) completed. "
            "See task record for summary and artifacts."
        )
