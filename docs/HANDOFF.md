# Handoff (2026-09-29)

State at the end of the session that built the demo, ran the pipelines on real data and
cut v0.0.1. Read `AGENTS.md`, `REPRODUCTION.md` and `TASKS.md` first; this file only adds
what they do not say.

## Where things are

- `main` is green (`make check`). No open PRs. Latest release: v0.0.1 (agents release only
  `0.0.x`; the old v0.1.0 tag and release were deleted on owner request).
- Real-data outputs and the exact commands are in `task-runs/20260929-real-census/`
  (`receipt.md` has the index and the numbers).
- Demo: `make demo` writes `demo-run/replay.html` (2D replay) and a text walkthrough.
  The README GIF is a recording of that page; capture mode is `replay.html#step=N&p=0..1`.

## Box state (`agent-dev-01`, user `robot`, shared disk 309 GB, usually 91-100% full)

Not in git, safe to delete once the receipt is no longer being re-counted:

- `~/runs/real/results/rnaseq/bam/SRR19152{324,325,327}.sorted.bam` (about 2.4 GB), the
  Bowtie2 index next to them, and `~/runs/real/data/` (RefSeq viral GenBank flat file, SA1 and
  host references, Pfam profiles).
- `~/runs/art-harness-20260929`: a clone used to run the wrappers.
- Tools: micromamba env `art` (`export MAMBA_ROOT_PREFIX=~/tools/mamba; eval "$(~/tools/mm/bin/micromamba shell hook -s bash)"; micromamba activate art`).

## Next actions, in order

1. **RNA-seq array share (GAP-11).** 4.5% and 3.8% at 15 min against the paper's 8%. First
   test: assign each fragment by its midpoint (paper-notes.md l.171-173) instead of
   featureCounts overlap; the 3 BAMs on the box are enough for that. Then run the other 9
   libraries (SRR19152326 failed on a full disk) on a dedicated volume.
2. **Live pilot** (owner decision needed: budget, model = Fable 5 in place of Mythos 5):
   1-2 tasks through the stage chain plus the L1 benchmark with a token cap.
3. **myRT profiles (GAP-10).** Census used 6 public Pfam RT profiles, not the 45 myRT
   profiles; AbiK is therefore not retained.
4. Rename the local checkout directory `art-harness` to `art-discovery-repro` from a fresh
   session (nothing in the repo depends on the name).
5. Owner closeout report.

## Gotchas learned

- A hook blocks git in the main checkout: work in `.worktrees/<task>-<YYYYMMDD>`. A fresh
  worktree has no dev deps, so run `make setup` there before the first commit or the
  commit hook reports "Tests failing".
- Run `git push -q -u origin HEAD` as its own command (the hook rejects it inside a
  compound command).
- On the box use `pkill -x <name>`, never `pkill -f` (it kills your own ssh session).
- The shared disk fills up during RNA-seq: keep `KEEP_INTERMEDIATES=0`, the SRA route
  (`FASTQ_SOURCE=sra`; ENA is about 100x slower from this box) and check `df` before each
  library. Other lanes also write to that disk.
- Local `main` in the main checkout can lag `origin/main`; fetch and compare before trusting it.
