# ART harness demo: what happened, and why

Agents: **rule policy (no model)**. Wall time 157s. Everything below was read back from `campaign/` after the run.

## What this is
A miniature of the campaign in Yoon et al. 2026 (Anthropic). The paper had ~950 model sessions search 1.9B protein clusters and find ART, a tandem-repeat-array reverse transcriptase system. This run keeps the paper's **workflow** and shrinks the **world** to nine synthetic contigs with a hidden, planted answer, and the **model** to a rule policy or a 0.6B local model.

## Why the workflow is shaped this way
| Harness element | What it does | Why (paper) | Seen in this run |
| --- | --- | --- | --- |
| Stage chain + scripted gates | five stages; a stage's tasks open only after the previous gate passes | model output is not trusted to say "done" (p.28) | 5 gates, all checked by code |
| Worker | plans, runs tools, writes summary | the agent that does the science | 10 sessions |
| Supervisor | accepts or returns for revision after checking claims against artifacts | independent review; 49 of 119 paper tasks were revised | 10 sessions, 1 task(s) revised |
| Curator + knowledge base | writes findings into a store later agents read | agents do not repeat or contradict earlier work | 9 entries |
| Triage | follow-ups are released or rejected with a written reason | keeps agent-proposed work bounded | 4 released, 1 rejected |
| Editor | reviews the report before it is filed | unsupported claims are caught before they leave | report filed |
| Git task records | every state change is a commit | the audit trail humans monitor | see the log below |

## The five stages
| stage | why it exists | tasks | status |
| --- | --- | --- | --- |
| 1_input_assembly | know exactly what is being analysed before anything else | t0001 | curated |
| 2_database_sweep | cast the widest net first (paper: 52 RT profile HMMs over 1.9B clusters) | t0002 | curated |
| 3_rt_classification | sort candidates so later stages compare like with like (paper: 9 classes) | t0003 | curated |
| 4_neighborhood_census | look next to every RT for an array (paper: 10,983 loci) | t0004, t0005, t0006, t0007, t0008 | curated |
| 5_deep_dives | spend effort only on what the census flagged (paper: 16 deep dives) | t0009 | completed |

## Every task
| task | stage | origin | revisions | status | tool headline |
| --- | --- | --- | --- | --- | --- |
| t0001 | 1_input_assembly | seed | 0 | curated | `contigs=9 total_nt=45440` |
| t0002 | 2_database_sweep | seed | 0 | curated | `rt_orfs=9 contigs_with_rt=9` |
| t0003 | 3_rt_classification | seed | 0 | curated | `classified=9 fragments=1 classes=A:2,B:4,C:2,fragment:1` |
| t0004 | 4_neighborhood_census | seed | 1 | curated | `assessed=8 arrays=4 no_array=3 not_assessed=1` |
| t0005 | 4_neighborhood_census | follow_up of t0004 | 0 | curated | `deep_dives=1 delimited=1 loci=L01` |
| t0006 | 4_neighborhood_census | follow_up of t0004 | 0 | curated | `deep_dives=1 delimited=1 loci=L02` |
| t0007 | 4_neighborhood_census | follow_up of t0004 | 0 | curated | `deep_dives=1 delimited=1 loci=L03` |
| t0008 | 4_neighborhood_census | follow_up of t0004 | 0 | curated | `deep_dives=1 delimited=1 loci=L09` |
| t0009 | 5_deep_dives | seed | 0 | completed | `deep_dives=0 delimited=0 loci=-` |

## A revision, in full
Task t0004 was returned 1 time(s). The supervisor's first reason (`t0004/verdict-r1.md`):

```
VERDICT: revise (revise)

- no tool call succeeded in this task, so no claim in the summary has evidence; run the tools named in the brief
```

## Triage rejections (each carries a written reason)
```
REJECTED AT TRIAGE
reason: arrays.tsv calls L06 'not_assessed'; a deep dive needs an array call

brief:
Deep dive locus L06: the scan could not assess it, so look for an array anyway.
```

## Knowledge base
- `t0001-findings-from-t0001-1-input-assembly.md`: Findings (rule agent copies the verified headlines): - contigs=9 total_nt=45440
- `t0002-findings-from-t0002-2-database-sweep.md`: Findings (rule agent copies the verified headlines): - rt_orfs=9 contigs_with_rt=9
- `t0003-findings-from-t0003-3-rt-classification.md`: Findings (rule agent copies the verified headlines): - classified=9 fragments=1 classes=A:2,B:4,C:2,fragment:1
- `t0004-findings-from-t0004-4-neighborhood-census.md`: Findings (rule agent copies the verified headlines): - assessed=8 arrays=4 no_array=3 not_assessed=1
- `t0005-findings-from-t0005-followup-by-worker.md`: Findings (rule agent copies the verified headlines): - deep_dives=1 delimited=1 loci=L01
- `t0006-findings-from-t0006-followup-by-worker.md`: Findings (rule agent copies the verified headlines): - deep_dives=1 delimited=1 loci=L02
- `t0007-findings-from-t0007-followup-by-worker.md`: Findings (rule agent copies the verified headlines): - deep_dives=1 delimited=1 loci=L03
- `t0008-findings-from-t0008-followup-by-worker.md`: Findings (rule agent copies the verified headlines): - deep_dives=1 delimited=1 loci=L09
- `t0009-findings-from-t0009-5-deep-dives.md`: Findings (rule agent copies the verified headlines): - deep_dives=0 delimited=0 loci=-

## Accounting next to the paper
| | this run | paper |
| --- | --- | --- |
| tasks | 9 | 119 |
| sessions | 31 | 949 |
| tasks revised >= 1x | 1 | 49 |
| triage rejections | 1 | 10 |
| reports filed | 1 | 19 |
| model tokens (in/out) | 0/0 | 11.3M/14.9M (+189.5M cache writes) |

## Scored against the hidden truth
- planted arrays found: **4/4**
- false array calls: **0** 
- decoy `phage_04 no_array`: handled correctly
- decoy `phage_05 no_array`: handled correctly
- decoy `phage_06 short_upstream`: handled correctly
- decoy `phage_07 dense_decoy`: handled correctly
- decoy `phage_08 fragment`: handled correctly

## The filed report
```markdown
# Array-associated RT candidates in the toy world

Demo campaign; the brief is a reconstruction, not Anthropic's text. Tasks: t0001, t0002, t0003, t0004, t0005, t0006, t0007, t0008, t0009.

## Result
4 array-positive loci among 8 assessed full-length RTs (9 RT candidates in total).

| locus | contig | class | R | deep dive | conserved block nt | median spacing nt |
| --- | --- | --- | --- | --- | --- | --- |
| L01 | phage_01 | B | 7 | delimited | 10 | 152 |
| L02 | phage_02 | B | 5 | delimited | 10 | 167 |
| L03 | phage_03 | C | 9 | delimited | 10 | 138 |
| L09 | phage_09 | A | 4 | delimited | 10 | 179 |

## Set aside
- fragments (under 250 aa): L08
- not assessed (too little upstream sequence): L06
- assessed, no array: L04, L05, L07

## Method
Motif ORF search, then the paper's k-mer array scan (Methods p.32 step 1, 100 mononucleotide shuffles) and array delimitation (step 2).

## Limits
Nine synthetic contigs and a motif in place of the 52 RT profile HMMs. No experiment here supports any biological claim.
```

## Audit trail (git log of the task records, oldest first)
```
task(t0001): open  [seed]
task(t0001): worker submitted summary
task(t0001): supervisor accepted
task(t0001): curated
task(t0002): open  [seed]
task(t0002): worker submitted summary
task(t0002): supervisor accepted
task(t0002): curated
task(t0003): open  [seed]
task(t0003): worker submitted summary
task(t0003): supervisor accepted
task(t0003): curated
task(t0004): open  [seed]
task(t0004): worker submitted summary
task(t0004): supervisor revise (1/3)
task(t0005): open followup by worker [follow_up]
task(t0006): open followup by worker [follow_up]
task(t0007): open followup by worker [follow_up]
task(t0008): open followup by worker [follow_up]
task(t0004): worker submitted summary
triage(t0004): reject follow-up — arrays.tsv calls L06 'not_assessed'; a deep dive needs an array call
task(t0004): supervisor accepted
task(t0004): curated
task(t0005): worker submitted summary
task(t0005): supervisor accepted
task(t0005): curated
task(t0006): worker submitted summary
task(t0006): supervisor accepted
task(t0006): curated
task(t0007): worker submitted summary
task(t0007): supervisor accepted
task(t0007): curated
task(t0008): worker submitted summary
task(t0008): supervisor accepted
task(t0008): curated
task(t0009): open  [seed]
task(t0009): worker submitted summary
task(t0009): supervisor accepted
task(t0009): curated
report(t0009): filed after editor review
task(t0009): report filed
```

## What is real and what is a stand-in
| real | stand-in |
| --- | --- |
| the orchestrator, gates, triage, roles, records, ledger | 9 synthetic contigs |
| the paper's k-mer array scan and array delimitation (`artharness.arrays`) | an RT motif instead of the 52 profile HMMs |
| the supervisor's claim check (code) | rule policy (no model) instead of Claude Mythos 5 |

## Try it
```
python -m artharness.demo                      # rule policy, offline
python -m artharness.demo --agent ollama      # small local model (qwen3:0.6b)
python -m artharness.demo --no-stage-revision # skip the staged sloppy pass
```
