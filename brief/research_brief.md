# Research brief (RECONSTRUCTION — not the verbatim original)

> **Provenance warning.** Anthropic's verbatim research brief (Supplementary Note 1)
> has not been published (see REPRODUCTION.md GAP-1). This document reconstructs the
> brief's skeleton from six anchors the paper states explicitly — five quoted, and
> anchor 5 cited as a NEGATIVE anchor (the verbatim brief's silence on repeats and
> arrays is itself the load-bearing fact). It is written to be *sufficient to run
> this harness*, not to impersonate Anthropic's text. Any campaign receipt must
> record: `brief: reconstruction-v0 (6 anchors), not verbatim`.
> Anchors: `docs/paper-notes.md` §External anchor points.

---

# Research brief: discovery of novel reverse-transcriptase systems

## Mission

Survey public metagenomic and genomic sequence data to **identify novel RT systems on
the basis of new partner-gene associations** [anchor 1]. The mission is directed at
**partner protein-coding genes**; non-coding DNA features are out of scope for
candidate promotion [anchor 2].

Background: reverse transcriptases (RTs) copy RNA into DNA. Bacterial and phage RTs
act inside larger systems, almost always with dedicated partner components — e.g.
retrons (RT + ncRNA + effector), diversity-generating retroelements (RT + accessory
components), abortive-infection RTs, and RT–CRISPR fusions. New
partner associations are how new RT systems have historically been recognized.

Reference systems to calibrate against (named examples):
**Ec86 retron (M24408.1), LtrA (group II intron, U01385), BPP-1 Brt (DGR),
AbiK (abortive infection), and an RT-Cas1 fusion** [anchor 4].

## Stages

The campaign proceeds in five stages; each stage closes with a scripted completion
check, and no later-stage task opens before the prior stage closes [anchor 3]:

1. **Input assembly** — assemble the working sequence database and its programmatic
   interface (predicted proteins, gene neighborhoods, biosample metadata,
   precomputed domain annotations; homology-search access to a reference protein set).
2. **Database sweep** — identify all RT candidates in the database using profile
   search against RT family profiles.
3. **RT classification** — classify candidates into known RT classes with a
   tiered scheme; classify thresholds against the reference systems above
   [anchor 6: thresholds "fixed beforehand" on the references + labeled RT sets].
4. **Neighborhood census** — sample loci spanning every RT class; score recurring
   protein families in RT neighborhoods as candidate partners against
   pre-committed selection criteria and pre-designated controls [anchor 6: caps and
   filters are yours to define before the census runs].
5. **Deep dives** — for each promoted candidate family, re-apply the promotion
   criteria, assemble alignments, and search annotations and the literature for
   prior descriptions. A family is confirmed when the evidence upholds it as a new
   RT partner; otherwise reject it or set it aside with reasons.

## Standing instructions

- Work from primary data. Where an observation is quantitative (copy counts,
  distances, score distributions), recount it from retrieved sequence.
- Test novelty claims by literature search with precise descriptors before reporting.
- If observations suggest a question that deserves a dedicated task, propose a
  follow-up task with a one-paragraph brief; do not silently expand scope.
  Follow-ups are triaged and may be rejected with a written reason.
- Record every finding in the shared knowledge base with task provenance.
- File candidate systems as reports: state the claim, the evidence, the confidence,
  and the experiment that would test the main prediction.

---

## GAP note (what the real brief may contain that this reconstruction does not)

Tone, hedges, preferred tool references, per-stage completion definitions, the exact
list of named reference accessions, any explicit anti-contamination instructions, and
the benchmark-generation instructions. Any campaign receipt must state:
`brief: reconstruction-v0 (6 anchors), not verbatim`.
