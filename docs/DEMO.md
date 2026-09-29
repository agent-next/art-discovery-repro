# Running demo: the whole workflow at toy scale

`python -m artharness.demo` runs the paper's campaign workflow end to end in about a
minute, offline, and writes `WALKTHROUGH.md` explaining what happened from its own
records. `--agent ollama` swaps the rule policy for a small local language model
(default `qwen3:0.6b`) that plays the same roles through the same prompts.

```bash
make demo                                   # rule policy, no model, deterministic
make demo-model                             # small model via Ollama (ollama serve; ollama pull qwen3:0.6b)
python -m artharness.demo --out /tmp/run --force --concurrency 4
```

Output lands in `demo-run/`: live commentary on the console, then
`demo-run/WALKTHROUGH.md`, `campaign/` (git-versioned task records, knowledge base,
ledger, shared tables) and `world/` (the input FASTA and the hidden truth).

## What it is

| | Paper (Yoon et al. 2026) | This demo |
| --- | --- | --- |
| World | 1.9B protein clusters, 10,983 RT loci | 9 synthetic contigs with a planted answer |
| Agents | Claude Mythos 5 in Claude Code, ~950 sessions | a rule policy, or a 0.6B local model |
| Stages | 5, closed by scripted gates | the same 5, same gate semantics |
| Tasks | 119, 49 revised, 10 triage rejections | ~9, 1 staged revision, 1 rejection (rules) |
| Science | HMM census, k-mer array scan, delimitation | motif ORF search, then the paper's own `kmer_scan` and `delimit_array` |
| Ending | 19 reports, tournament | one report, editor review, scored against the truth |

The orchestrator, gates, triage, roles, record store, knowledge base and ledger are the
production code in `src/artharness/`, not a copy. Only the world, the two agents and the
five tools live in `src/artharness/demo/`.

## The world

Nine contigs, each with one candidate RT gene (motif `[YF][ADV]DD`, NOT-IN-PAPER: the
paper uses 52 profile HMMs). Planted: four true tandem-repeat arrays upstream of an RT,
and four traps a careless analysis would fall into: an RT fragment under 250 aa, an RT
with too little upstream sequence to assess, an RT with a dense tandem repeat (spacing
under 100 nt, which the paper's suppression rule must reject), and plain RTs with no
array. The truth is written to `world/truth.json`; the tools never read it, only the
scorecard does.

## How one task flows

```
seed brief / released follow-up
        v
   worker: PLAN -> TOOL lines -> harness runs tools (tool_log.jsonl) -> SUMMARY + HEADLINE lines
        v
   completion check: summary.md exists (else worker re-dispatched, stall after 3 here, 10 in the paper)
        v
   triage: worker/supervisor follow-ups released or rejected with a written reason
        v
   supervisor: scripted claim check + verdict.  revise -> verdict-rN.md -> worker sees it -> repeat
        v
   curator: findings -> knowledge base (later prompts include relevant entries)
        v
   stage gate (after every task of the stage is drained): code, not opinion
```

## Why each piece is there

- **Scripted gates.** A model saying "done" is not evidence. Each gate recomputes what it
  needs from the raw input (`campaign.make_gates`). Kill a stage's output and the chain
  stops, as in the small-model run below.
- **Claim check.** The supervisor's duty in the paper is to check quantitative claims
  against artifacts. Here `agents.evidence_problems` requires each tool `HEADLINE:` line
  quoted verbatim and rejects any number of two or more digits that no tool printed. For
  the small model this check is a necessary condition of "accept": it can be talked into
  nothing.
- **Scope.** A worker may only call the tools its brief names; anything else must become
  a follow-up, which is what the paper's worker prompt says ("do not expand the current
  task's scope"). A tool call outside the brief returns an error and is logged.
- **Triage.** Follow-ups are cheap for a model to invent. Each is released or rejected
  with a reason written into the task record.
- **Revision feedback.** A revision pass receives the supervisor's last verdict
  (`Roles.worker`), so "revise" carries information.
- **Git records.** Every state change is a commit; the audit trail a human monitors is
  `git log` in `campaign/`.

## What the small model does (qwen3:0.6b, recorded run)

A 0.6B model is a poor scientist and a good demonstration of why the harness has its
structure. In the recorded run (`task-runs/20260929-demo/console-ollama.txt`, 6 tasks,
about 2 minutes on CPU) it:

- called the right tool every time, but 5 of 6 tasks were returned for revision: it
  dropped or reworded the `HEADLINE:` line, or stated "100" (a number no tool printed).
  The claim check named the exact missing line and every task passed on a later attempt;
- copied the protocol's made-up example (`HEADLINE: widgets=3 sprockets=1`) as a
  follow-up, and once wrote `<no follow-up>`; triage rejected both with a reason (4
  rejections in all);
- proposed 1 of the 4 possible deep dives. The stage-5 task ran the rest, and the stage-5
  gate checked that every array had one;
- finished with the report filed and 4/4 planted arrays found, 0 false calls.

In earlier development runs the same model also tried a tool outside its brief (the
scope guard refused) and, with a looser protocol, stalled 2 tasks and left a gate
failing, which stopped the campaign before the report: the gates do stop it. Runs vary
with seed and prompt. Nothing is scripted in the model path; the tests drive it with an
HTTP stub that speaks the Ollama API.

## Honest limits

- It demonstrates the workflow, not the science: nine synthetic contigs prove nothing
  about biology. The deep-dive "conserved block" is the 10-nt seed the paper's
  delimitation reports, not the full repeat.
- The rule policy's stage-4 first pass is deliberately sloppy (labelled in the task
  record) so the revise loop appears; `--no-stage-revision` removes it.
- Token counts for the rule policy are 0 (no model). The model path reports Ollama's
  prompt and output token counts; the paper's cache-write class does not exist here.
- No launch session is run (the paper had one); the caller supplies the briefs.
- A non-local model server needs `ARTHARNESS_ALLOW_LIVE=1`, the same owner gate as the
  Claude backend.
