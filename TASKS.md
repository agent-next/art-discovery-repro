# TASKS.md — art-harness build queue

Goal: max reproduction of Yoon et al. 2026 (ART harness). Push every layer until a
real blocker is hit; every blocker has a WHY/HOW/WHAT entry in REPRODUCTION.md §5.

## Done
- [x] Repo scaffold, packaging, CI; `make check` gate (ruff + shellcheck + pytest with coverage floor)
- [x] Harness core: records / knowledge / roles / orchestrator / tournament / accounting / backends
- [x] `arrays.py`: k-mer scan, delimitation, PWM extension, cross-scan (Methods p.32)
- [x] Pipeline wrappers (census 01–06, art_family, rnaseq, db subset), connectors, 12 skill guides, benchmark
- [x] `docs/paper-notes.md`, `REPRODUCTION.md`, `brief/research_brief.md` (6-anchor reconstruction)
- [x] Five-round independent review; v0.1.0 tagged; repo renamed `art-discovery-repro`
- [x] Stage gates evaluate after each stage's work; follow-up dedupe; eval removed from pipeline wrappers

## Open
- [ ] Real-data runs of the pipeline wrappers: install the bio tools (diamond, mmseqs2,
      hmmer, mafft, iqtree, bowtie2, seqkit) in a scratch env, then run MarsHill scan,
      census subset, RNA-seq reanalysis (PRJNA836150) and record receipts under `task-runs/`
- [ ] Live pilot: 1–2 tasks through the stage chain + L1 benchmark with a token cap
      (owner decision: budget, model = Fable 5 in place of Mythos 5)
- [ ] Threaded dispatch up to `max_concurrent_sessions` (paper: ≤58) — currently sequential
- [ ] Contig-length-aware cluster representative in `family_definition.sh` (GAP in script header)
- [ ] Wire the GenBank BLASTP branch into the ART pool (output currently unused)
- [ ] Rename the local checkout directory `art-harness` → `art-discovery-repro` from a fresh session
- [ ] Owner closeout report (goal-html-report)

## Blockers ledger
See REPRODUCTION.md §5 (GAP-1 … GAP-7): unreleased brief and transcripts, Mythos 5
internals, planetary database scale, wet lab, replicate-campaign budget.
