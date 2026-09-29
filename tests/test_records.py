from pathlib import Path

from artharness.knowledge import KnowledgeBase
from artharness.records import RecordStore, TaskOrigin, TaskStatus


def test_task_lifecycle_roundtrip(tmp_path: Path):
    store = RecordStore(tmp_path, use_git=False)
    rec = store.create("t0001", "find novel RT partners", "4_neighborhood_census",
                       TaskOrigin.SEED, label="census")
    assert rec.status is TaskStatus.OPEN
    rec.status = TaskStatus.EXECUTED
    rec.revisions = 2
    store.update(rec, "task(t0001): executed")
    got = store.get("t0001")
    assert got.status is TaskStatus.EXECUTED and got.revisions == 2
    assert got.origin is TaskOrigin.SEED
    store.write_text("t0001", "summary.md", "found 3 partner families")
    assert "partner families" in store.read_text("t0001", "summary.md")
    assert store.artifact_dir("t0001").is_dir()
    assert [r.task_id for r in store.list_tasks()] == ["t0001"]
    assert store.next_task_id(7) == "t0007"


def test_create_rejects_duplicate(tmp_path: Path):
    store = RecordStore(tmp_path, use_git=False)
    store.create("t0001", "b", "5_deep_dives", TaskOrigin.DEEP_DIVE)
    import pytest

    with pytest.raises(ValueError):
        store.create("t0001", "b", "5_deep_dives", TaskOrigin.DEEP_DIVE)


def test_knowledge_retrieval_and_context_block(tmp_path: Path):
    kb = KnowledgeBase(tmp_path / "kb")
    kb.add("t0001", "retron RTs act in phage defense",
           "Retron RTs commonly act in phage defense systems with an ncRNA upstream.")
    kb.add("t0002", "CRISPR arrays make crRNAs",
           "CRISPR arrays are transcribed into crRNA guides from repeat-spacer units.")
    hits = kb.relevant("does the retron have an ncRNA upstream?", limit=1)
    assert len(hits) == 1 and "t0001" in hits[0][0].name
    block = kb.context_block("crRNA guides from arrays")
    assert "CRISPR arrays" in block and "relevance" in block
    assert kb.context_block("completely unrelated quantum coffee") == ""


def test_git_backed_store_commits_the_trail(tmp_path):
    # Audit finding 5: the default use_git=True path had ZERO tests — every
    # test in the suite disabled it, so a broken commit() would ship green.
    # Real oracle: a real git repo (init in tmp), a task transition, then
    # `git log`/`git show` must contain the transition in a commit.
    import subprocess

    root = tmp_path / "camp"
    root.mkdir()
    # -b campaign: the host's global git-guard refuses commits on a branch
    # named 'main' (fresh inits default to it); the store's trail mechanics
    # are the test subject, not the guard
    subprocess.run(["git", "init", "-q", "-b", "campaign"], cwd=root, check=True)
    subprocess.run(["git", "config", "user.email", "t@example.com"], cwd=root, check=True)
    subprocess.run(["git", "config", "user.name", "t"], cwd=root, check=True)

    store = RecordStore(root, use_git=True)
    assert store.use_git is True  # inside a repo -> git trail ON (the default)
    rec = store.create("t0001", "brief", "1_input_assembly", TaskOrigin.SEED)
    store.update(rec, "completed after gate")
    out = subprocess.run(["git", "log", "--oneline"], cwd=root,
                         capture_output=True, text=True, check=True)
    assert out.stdout.strip(), "no commits recorded"
    show = subprocess.run(["git", "show", "--stat", "HEAD"], cwd=root,
                          capture_output=True, text=True, check=True)
    assert "records/" in show.stdout  # the trail file itself is committed
    # outside a repo the store degrades, it does not die (self.use_git False)
    bare = RecordStore(tmp_path / "norepo", use_git=True)
    assert bare.use_git is False
