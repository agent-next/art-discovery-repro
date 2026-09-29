"""A revision pass must show the worker why the supervisor returned the last one
(paper p.28: supervisors "returned them for revision")."""

from pathlib import Path

from artharness.knowledge import KnowledgeBase
from artharness.records import RecordStore, TaskOrigin
from artharness.roles import Roles
from artharness.runner.base import ScriptedBackend


def _roles(tmp_path: Path):
    store = RecordStore(tmp_path, use_git=False)
    rec = store.create("t0001", "count things", "1_input_assembly", TaskOrigin.SEED)
    backend = ScriptedBackend()
    return Roles(backend, store, KnowledgeBase(tmp_path / "kb"), None), store, rec, backend


def test_first_pass_has_no_feedback_section(tmp_path: Path):
    roles, _, rec, backend = _roles(tmp_path)
    roles.worker(rec)
    assert "supervisor feedback" not in backend.calls[-1].user_prompt


def test_revision_pass_carries_the_latest_verdict_only(tmp_path: Path):
    roles, store, rec, backend = _roles(tmp_path)
    store.write_text("t0001", "verdict-r1.md", "VERDICT: revise\n\nfirst complaint\n")
    store.write_text("t0001", "verdict-r2.md", "VERDICT: revise\n\nsecond complaint\n")
    store.write_text("t0001", "verdict-r10.md", "VERDICT: revise\n\ntenth complaint\n")
    rec.revisions = 10
    roles.worker(rec)
    prompt = backend.calls[-1].user_prompt
    assert "tenth complaint" in prompt
    assert "first complaint" not in prompt and "second complaint" not in prompt
