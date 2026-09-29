# art-harness

Faithful, self-contained reproduction of the autonomous research harness from:

> Yoon, Athukoralage, Ameisen, Kauderer-Abrams, Perry, Durrant et al. (Anthropic),
> **"Autonomous AI agents discover reverse transcriptases with tandem repeat arrays"** (2026-09-23).
> Paper PDF: https://www-cdn.anthropic.com/22573675ada52a8ca8a97a1a4b4326b2f208a071.pdf
> Announcement: https://www.anthropic.com/news/claude-discovers-novel-enzyme-system

In that paper, ~950 Claude Code sessions (119 tasks, 215.6M tokens, 21.5 h) surveyed
reverse-transcriptase loci across 1.9B protein clusters and discovered **ART**
(array-associated reverse transcriptases), a novel CRISPR-like system in jumbo phages.

This repo rebuilds the three layers of that system with public tools:

| Layer | Paper component | Here |
| --- | --- | --- |
| **Harness** | launch/worker/supervisor/curator/editor sessions, task records in version control, scripted stage gates, triage with written rejections, shared knowledge base, token accounting | `src/artharness/` |
| **Science pipeline** | RT census (52 HMMs → 198,290 clusters → 9 classes → 10,983 loci → 16 deep dives) + ART family/array/phylogeny/RNA-seq analyses, all paper-exact parameters | `pipeline/` |
| **Benchmark** | fixed-input L1–L5 benchmark (3,500 attempts in the paper), 10-claim rubric, offline grader (paper rule; the paper's LLM judge is NOT-IN-PAPER as a mechanism here — benchmark/rubric.py), report tournament (Bradley–Terry) | `benchmark/` |

## Status

v0.1.0 tagged. All three layers are implemented and tested offline (`make check`:
ruff, shellcheck, pytest with a coverage floor). No live model campaign, benchmark
run or real-data pipeline run has been executed yet: those are budget/compute-gated.
`REPRODUCTION.md` maps every paper element to its component and state and lists the
unreproducible gaps; `TASKS.md` holds the open work.

## Safety gates

- The LLM session backend is pluggable (`src/artharness/runner/`). Tests run fully offline
  against a scripted backend. **Launching a real LLM campaign is an owner decision** — the
  harness never calls a paid model API unless a backend is explicitly configured.
- All wet-lab steps from the paper are documented in `pipeline/` as protocols only;
  nothing here operates lab equipment.

## Check

```bash
make check   # ruff + pytest
```
