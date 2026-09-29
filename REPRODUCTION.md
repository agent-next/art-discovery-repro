# REPRODUCTION.md — paper element → component map, and the gaps ledger

Paper: Yoon et al. (Anthropic) 2026, "Autonomous AI agents discover reverse
transcriptases with tandem repeat arrays". Facts ledger: `docs/paper-notes.md`.

## 1. Map

| # | Paper element | Paper ref | Component | State |
|---|---|---|---|---|
| 1 | Launch agent, 5-stage chain, scripted completion gates | Methods p.28 | `orchestrator.run_stage_chain` + caller-supplied `Gate` callables | done |
| 2 | Task records in version control, readable by every agent | Methods p.28 | `records.RecordStore` (git-backed tree; git optional for tests) | done |
| 3 | Worker: plan → execute → summary + artifacts; proposes follow-ups | Methods p.28 | `roles.Roles.worker` + role prompts | done (prompt text is ours; original not public — GAP-1) |
| 4 | Supervisor accept/revise review; writes follow-up briefs | Methods p.28, p.31 | `Roles.supervisor`; orchestrator revision loop; `verdict-rN.md` | done |
| 5 | Stall rules: 10 revisions / 10 failed completion checks | Methods p.30 | `config.max_revisions/max_gate_failures`, `Orchestrator._dispatch` | done |
| 6 | Triage queue, release or reject with written reason | Methods p.28 | `Orchestrator.propose_followup` + `TriagePolicy`, `triage-rejection-*.md` | done |
| 7 | Curator → shared knowledge base → relevant entries into later prompts | Methods p.28 | `knowledge.KnowledgeBase` (term-overlap retrieval) | done (retrieval scheme NOT-IN-PAPER) |
| 8 | Editor review before filing reports | Methods p.28 | `Orchestrator.file_report`, `Roles.editor` | done |
| 9 | Concurrency ≤58, sandbox 60 CPU/192 GiB/no GPU | Methods p.28 | `config.CampaignConfig` | done (values pinned by test); dispatch is sequential — concurrency UNIMPLEMENTED |
| 10 | Connectors: protein DB, literature, KB, GPU queue | Methods p.28 | `connectors/` | done (offline-first clients; live lanes per connector header) |
| 11 | ~140 skills; 7 survey guides named | Methods p.28 | `skills/` (12 guides; full 140-library out of scope) | done (GAP-2) |
| 12 | Campaign accounting: sessions/roles, agent-hours, token classes | Methods p.30 | `accounting.SessionLedger` | done |
| 13 | Report tournament: 342 games, rubric weights, soundness auto-lose, BTL | Methods p.30 | `tournament.py` | done |
| 14 | RT census: 52 HMMs → 203,381 → 198,290 clusters → 9 classes | Methods p.29 | `pipeline/census/01–04` | done |
| 15 | Neighborhood sampling: 7,308 anchors → 10,983 loci | Methods p.29-30 | `pipeline/census/05` | done + tests (`test_pipeline_logic`) |
| 16 | Partner scoring: 3 filters, controls, 3,564 → 16 | Methods p.30 | `pipeline/census/06` | done + tests (`test_pipeline_logic`) |
| 17 | ART family definition: QQM14740.1 → 95 members | Methods p.31-32 | `pipeline/art_family/family_definition.sh` | done as a flag-exact wrapper (dry-run + flag tests); never run on real data; two GAPs noted in the script header |
| 18 | k-mer array scan (20×14-mer seeds, shuffles, R≥3) | Methods p.32 | `artharness.arrays.kmer_scan` | done + tests |
| 19 | Array delimitation (10-mer, 30% tol, 200/2000 shuffles, PWM) | Methods p.32 | `artharness.arrays.delimit_array` + `pwm_extend` + `cross_scan` | done + tests |
| 20 | Phylogeny: 774 set, MAFFT L-INS-i, IQ-TREE Q.pfam+F+R6 | Methods p.33 | `pipeline/art_family/phylogeny.sh` | done as a flag-exact wrapper (dry-run + flag tests); never run on real data |
| 21 | RNA-seq reanalysis: PRJNA836150, Bowtie2, 8% at 15 min | Methods p.36-37 | `pipeline/rnaseq/sa1_infection.sh` | done as a flag-exact wrapper (dry-run, flag and injection tests); never run on real data |
| 22 | Wet-lab protocols (plasmids, small-RNA-seq) | Methods p.36-37 | documented only — no lab (GAP-6) | documented |
| 23 | Benchmark: L1–L5, 3,500 attempts, 10-claim rubric, judge | Methods p.38 | `benchmark/` | done offline (levels, rubric, grader, runner tested); no live run (GAP-7 budget gate) |
| 24 | Replicate campaigns (10×) + transcript forensics | Methods p.38 | rerun via orchestrator + `experiments/forensics.py` (identifier search, ≥200-nt DNA + repeat-remark parsing) | forensics done; reruns budget-gated (GAP-7) |
| 25 | Interpretability (Evo2/gLM2 profiles; Mythos 5 sparse signals) | Methods p.38-39 | out of repo scope — needs model internals (GAP-3) | not reproducible |
| 26 | Serendipity chain t0010 → t0062 | Results p.3-5 | emerges from harness if workers/supervisors behave similarly | not guaranteed (paper: 0/10 reruns) |

## 2. Verification ladder (what "reproduced" means here)

1. **Unit green**: every module tested offline (`make check`).
2. **Dry-run fidelity**: pipeline wrappers print paper-exact commands (`--dry-run`).
3. **Data-anchored**: array code recovers planted synthetic ART arrays; on real data,
   `scripts/scan_genome.py` on the downloaded MarsHill genome (MW248466.1, ENA) must
   yield its 5-copy shuffle-controlled array call (seed ATATGAATACGTAT; recorded in
   `task-runs/20260924-devin-build/real-genome-scan.txt`). [Corrected 2026-09-24 S1
   review: this line previously claimed a 14-copy call from a nonexistent script.]
4. **Benchmark parity**: L1–L3 rerun with Fable 5; compare per-level recognition rates
   against the paper's ranges (L1 ≥90%, L3 model-dependent 32–96%).
5. **Campaign scale**: full 119-task campaign (GAP-4 compute + GAP-1 brief fidelity).
6. **Discovery**: NOT a reproduction target — the paper itself shows 0/10 reruns find
   the array. We reproduce the *possibility space*, not the accident.

## 3. Budget to rerun (from paper numbers)

- Full campaign ≈ 215.6M tokens (components 11.3M in + 14.9M out + 189.5M
  cache-write sum to 215.7M at display precision — the PDF's own rounding; the
  stated total 215.6M is what the paper prints),
  949 sessions, 21.5 h wall clock at ≤58 concurrent. At current frontier Claude API
  pricing: low-thousands of USD per campaign (estimate, not quoted).
- Benchmark (3,500 isolated attempts, ≤1M output tokens each) is the expensive half if
  run at paper scale; scale via `--attempts`.

## 4. Skills library

The paper names 7 metagenomic survey guides + libraries for database use, literature
access, GPU dispatch, report editing, bioinformatic tools. We ship those 12 guides
(`skills/survey/`, `skills/tools/`) as one-page method docs; the remaining ~128
guides' titles/contents were never published (GAP-2).

## 5. Gaps ledger — WHY / HOW / WHAT (the "cannot push further" list)

- **GAP-1 Research brief verbatim (Supp. Note 1).**
  WHY: not in the 40-page PDF (which contains Supp Figs 1–7 but no Notes); not linked
  on the news page; absent from HN (752 comments scraped), bioRxiv/arXiv/Semantic
  Scholar/OpenAlex (paper not yet indexed), Wayback snapshot.
  HOW: reconstructed from 6 anchors in the paper (`brief/research_brief.md`), each
  anchor cited; harness runs accept any brief as input.
  WHAT WOULD UNLOCK: Anthropic publishing the supplement (or the bioRxiv version).

- **GAP-2 Full skill library (~140 docs) and session transcripts.**
  WHY: never released; transcripts are referenced (Supp. Fig. 1A/2) but not deposited.
  HOW: ship the 12 named/derivable guides; prompt text ours.
  WHAT WOULD UNLOCK: Anthropic transcript release.

- **GAP-3 Model internals.**
  WHY: repeat-signal interpretability requires the Mythos 5 checkpoint (restricted
  access); public Fable 5 weights are not open either.
  HOW: behavioral substitutes — Evo 2/gLM2 likelihood profiles are computable from
  public checkpoints (queued in experiments/); benchmark L1–L5 measures the behavior
  without internals.
  WHAT WOULD UNLOCK: open-weights release or Anthropic-hosted probing access.

- **GAP-4 Planetary database build.**
  WHY: 15.8B proteins clustered with DIAMOND needs multi-TB RAM/disk and days of
  cluster time; not available in this workspace (229 GB free, single host).
  HOW: `pipeline/db/build_subset.sh` reproduces the exact recipe on any assembly
  subset; document scale-out path (DIAMOND linear-time clustering is embarrassingly
  parallel over partitions; the paper partitioned the database in 64).
  WHAT WOULD UNLOCK: a cluster allocation (Jetstream2/AutoDL is not sized for this;
  would need ≥2 TB RAM or a carefully partitioned disk-backed run).

- **GAP-5 Delimitation tail (PWM copy-extension step).** RESOLVED 2026-09-24:
  `pwm_extend` (log-odds PWM; copies = matches above the max of 200 block-shuffled
  regions) and `cross_scan` (every array's PWM scanned against every other locus to
  group arrays sharing a repeat) implemented and tested
  (`tests/test_arrays_pwm.py`). Pseudocount value NOT-IN-PAPER (1e-3 default).

- **GAP-6 Wet lab.**
  WHY: no laboratory; protocols require BSL-1/2 work, Illumina run, cloning.
  HOW: full reagent/protocol notes in `docs/paper-notes.md` §RNA-seq; public-data
  reanalysis (PRJNA836150) covers the computational half.
  WHAT WOULD UNLOCK: a partner lab; the paper's own wet work was human-run.

- **GAP-7 Replicate-campaign statistics.**
  WHY: 10 full campaigns ≈ 2.2B tokens ≈ real money; owner-gated.
  HOW: orchestrator supports reruns; forensics = identifier search over task records
  (≥200-nt DNA string scan) as specified in the paper.
  WHAT WOULD UNLOCK: owner budget approval.

## 6. Model delta

Paper model: Claude Mythos 5 (restricted twin of Fable 5). This repo's default
backend model id: Fable 5. Any campaign receipt must record
`model: <id> (paper: Mythos 5)` and treat numeric parity as approximate.
