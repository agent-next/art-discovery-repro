# ART harness demo: what happened, and why

Agents: **small model qwen3:0.6b via Ollama**. Wall time 76s. Everything below was read back from `campaign/` after the run.

## What this is
A miniature of the campaign in Yoon et al. 2026 (Anthropic). The paper had ~950 model sessions search 1.9B protein clusters and find ART, a tandem-repeat-array reverse transcriptase system. This run keeps the paper's **workflow** and shrinks the **world** to nine synthetic contigs with a hidden, planted answer, and the **model** to a rule policy or a 0.6B local model.

## Why the workflow is shaped this way
| Harness element | What it does | Why (paper) | Seen in this run |
| --- | --- | --- | --- |
| Stage chain + scripted gates | five stages; a stage's tasks open only after the previous gate passes | model output is not trusted to say "done" (p.28) | 5 gates, all checked by code |
| Worker | plans, runs tools, writes summary | the agent that does the science | 12 sessions |
| Supervisor | accepts or returns for revision after checking claims against artifacts | independent review; 49 of 119 paper tasks were revised | 12 sessions, 5 task(s) revised |
| Curator + knowledge base | writes findings into a store later agents read | agents do not repeat or contradict earlier work | 6 entries |
| Triage | follow-ups are released or rejected with a written reason | keeps agent-proposed work bounded | 1 released, 4 rejected |
| Editor | reviews the report before it is filed | unsupported claims are caught before they leave | report filed |
| Git task records | every state change is a commit | the audit trail humans monitor | see the log below |

## The five stages
| stage | why it exists | tasks | status |
| --- | --- | --- | --- |
| 1_input_assembly | know exactly what is being analysed before anything else | t0001 | curated |
| 2_database_sweep | cast the widest net first (paper: 52 RT profile HMMs over 1.9B clusters) | t0002 | curated |
| 3_rt_classification | sort candidates so later stages compare like with like (paper: 9 classes) | t0003 | curated |
| 4_neighborhood_census | look next to every RT for an array (paper: 10,983 loci) | t0004, t0005 | curated |
| 5_deep_dives | spend effort only on what the census flagged (paper: 16 deep dives) | t0006 | completed |

## Every task
| task | stage | origin | revisions | status | tool headline |
| --- | --- | --- | --- | --- | --- |
| t0001 | 1_input_assembly | seed | 0 | curated | `contigs=9 total_nt=45440` |
| t0002 | 2_database_sweep | seed | 1 | curated | `rt_orfs=9 contigs_with_rt=9` |
| t0003 | 3_rt_classification | seed | 1 | curated | `classified=9 fragments=1 classes=A:2,B:4,C:2,fragment:1` |
| t0004 | 4_neighborhood_census | seed | 1 | curated | `assessed=8 arrays=4 no_array=3 not_assessed=1` |
| t0005 | 4_neighborhood_census | follow_up of t0004 | 1 | curated | `deep_dives=1 delimited=1 loci=L01` |
| t0006 | 5_deep_dives | seed | 2 | completed | `deep_dives=3 delimited=3 loci=L02,L03,L09` |

## A revision, in full
Task t0002 was returned 1 time(s). The supervisor's first reason (`t0002/verdict-r1.md`):

```
VERDICT: revise (revise)

- the summary must quote the tool result verbatim; missing line: HEADLINE: rt_orfs=9 contigs_with_rt=9
VERDICT: accept

The summary correctly quotes the tool result, and there are no concrete errors to revise. The task is complete. 

VERDICT: accept
```

## Triage rejections (each carries a written reason)
```
REJECTED AT TRIAGE
reason: the follow-up names no locus, so there is nothing concrete to run

brief:
HEADLINE: widgets=3 sprockets=1
```
```
REJECTED AT TRIAGE
reason: the follow-up names no locus, so there is nothing concrete to run

brief:
HEADLINE: widgets=3 sprockets=1
```
```
REJECTED AT TRIAGE
reason: the follow-up names no locus, so there is nothing concrete to run

brief:
<no follow-up>
```
```
REJECTED AT TRIAGE
reason: the follow-up names no locus, so there is nothing concrete to run

brief:
<no such line>
```

## Knowledge base
- `t0001-findings-from-t0001-1-input-assembly.md`: - Contigs: 9 total (45440 nt) with GC content ranging from 0.417 to 0.444.   - Sprocket: Found and set aside for widget assembly.   - Task IDs: t0001, t0002.
- `t0002-findings-from-t0002-2-database-sweep.md`: - Found 9 RT-orfs across all six reading frames.   - Set aside 1 RT fragment.   - Tools: find_rt_orfs.
- `t0003-findings-from-t0003-3-rt-classification.md`: - Task t0003: Found 9 RT candidates with 2, 4, and 2 fragments classified. Tools used: classify_rt.
- `t0004-findings-from-t0004-4-neighborhood-census.md`: - **Foundings**: 4 arrays were assessed, with 3 not assessed and 1 assessed.   - **Evidence**: 4 arrays were analyzed, with no array not assessed.   - **Task IDs**: t0004, t0004, t0004.
- `t0005-findings-from-t0005-followup-by-worker.md`: - Found deep_dives=1 delimited=1 loci=L01.   - Evidence: Arrays analyzed, single locus delimited.   - Task ID: t0005.
- `t0006-findings-from-t0006-5-deep-dives.md`: - Every array-positive locus has been delimit and analyzed.   - Tools: deep_dive.   - Task ID: t0006.

## Accounting next to the paper
| | this run | paper |
| --- | --- | --- |
| tasks | 6 | 119 |
| sessions | 31 | 949 |
| tasks revised >= 1x | 5 | 49 |
| triage rejections | 4 | 10 |
| reports filed | 1 | 19 |
| model tokens (in/out) | 31530/3784 | 11.3M/14.9M (+189.5M cache writes) |

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

Demo campaign; the brief is a reconstruction, not Anthropic's text. Tasks: t0001, t0002, t0003, t0004, t0005, t0006.

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
triage(t0001): reject follow-up — the follow-up names no locus, so there is nothing concrete to run
task(t0001): worker submitted summary
task(t0001): supervisor accepted
task(t0001): curated
task(t0002): open  [seed]
task(t0002): worker submitted summary
task(t0002): supervisor revise (1/3)
task(t0002): worker submitted summary
task(t0002): supervisor accepted
task(t0002): curated
task(t0003): open  [seed]
task(t0003): worker submitted summary
task(t0003): supervisor revise (1/3)
task(t0003): worker submitted summary
task(t0003): supervisor accepted
task(t0003): curated
task(t0004): open  [seed]
triage(t0004): reject follow-up — the follow-up names no locus, so there is nothing concrete to run
task(t0004): worker submitted summary
task(t0004): supervisor revise (1/3)
task(t0005): open followup by worker [follow_up]
task(t0004): worker submitted summary
task(t0004): supervisor accepted
task(t0004): curated
triage(t0005): reject follow-up — the follow-up names no locus, so there is nothing concrete to run
task(t0005): worker submitted summary
task(t0005): supervisor revise (1/3)
task(t0005): worker submitted summary
task(t0005): supervisor accepted
task(t0005): curated
task(t0006): open  [seed]
task(t0006): worker submitted summary
task(t0006): supervisor revise (1/3)
triage(t0006): reject follow-up — the follow-up names no locus, so there is nothing concrete to run
task(t0006): worker submitted summary
task(t0006): supervisor revise (2/3)
task(t0006): worker submitted summary
task(t0006): supervisor accepted
task(t0006): curated
report(t0006): filed after editor review
task(t0006): report filed
```

## What is real and what is a stand-in
| real | stand-in |
| --- | --- |
| the orchestrator, gates, triage, roles, records, ledger | 9 synthetic contigs |
| the paper's k-mer array scan and array delimitation (`artharness.arrays`) | an RT motif instead of the 52 profile HMMs |
| the supervisor's claim check (code) | small model qwen3:0.6b via Ollama instead of Claude Mythos 5 |

## Try it
```
python -m artharness.demo                      # rule policy, offline
python -m artharness.demo --agent ollama      # small local model (qwen3:0.6b)
python -m artharness.demo --no-stage-revision # skip the staged sloppy pass
```
