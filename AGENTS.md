# AGENTS.md — art-harness

GitHub repo: `agent-next/art-discovery-repro` (renamed from `art-harness`; the local
checkout directory keeps the old name).

Org rules: https://github.com/agent-next/.github/blob/main/AGENT-STANDARD.md (hard limits, PR/merge policy).

## What this repo is

A reproduction of the Anthropic ART autonomous-research-harness paper (see README).
Parameter sources are cited as `paper p.NN` / `paper Methods §name` — every constant in
`pipeline/` and `benchmark/` traces to the paper PDF (local copy referenced in
`docs/paper-notes.md`).

## Orient

- `git status --short`, `git branch --show-current`, `git worktree list`,
  `gh pr list --state open`.
- Read first: `README.md`, `REPRODUCTION.md` (paper-element → component map),
  `TASKS.md` (build queue).

## Setup

`make setup` creates `.venv` and installs the package with dev extras
(`pip install -e '.[dev]'`) — the same install CI runs. Requires python3 >= 3.11
(CI pins 3.12). `uv.lock` is committed; `uv sync --extra dev` gives the pinned
equivalent.

## Check

`make check` — the gate `.github/workflows/ci.yml` runs: `ruff check src tests
pipeline benchmark` + `shellcheck -S warning pipeline/*/*.sh` + `python -m pytest`
with a coverage floor (offline; no network, GPU, or paid APIs). The Makefile uses
`.venv/bin/` automatically when it exists, so no `source .venv/bin/activate` is needed.
Single test file: `.venv/bin/python -m pytest tests/test_<name>.py`.

## Rules for agents working here

- Read `REPRODUCTION.md` before touching anything: it maps every paper element to code.
- Python under `src/artharness/` and `pipeline/*.py` is tested code: `make check` must
  pass before any commit. Shell scripts under `pipeline/*.sh` are parameter wrappers —
  keep every flag identical to the paper citation next to it.
- Do not invent paper parameters. If a value is not in the paper, mark it
  `# NOT-IN-PAPER: <reason>` in code and `GAP:` in docs.
- The reconstructed research brief lives in `brief/research_brief.md`; it is a
  reconstruction from six anchors (cited inline), NOT the verbatim original.
  Never present it as Anthropic's text.
- Task receipts: `task-runs/<YYYYMMDD>-<topic>/`.
- Commit small, one logical change per commit, Conventional Commits, no AI attribution.
- Worktrees for multi-file work: `.worktrees/<task>-<YYYYMMDD>` (org convention).
- Offline by default: unit tests must not hit the network. Live connectors get a
  `--live` flag and a recorded-fixture test mode.

## Done

Branch per change -> PR. New code gets a test with a real oracle. `make check`
green locally; CI green before merge. Receipts (commands + real output) go in
the PR body.
