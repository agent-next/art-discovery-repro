"""Repeat-array detection and delimitation — paper Methods p.32, parameter-exact.

Step 1 (k-mer scan, "Repeat-array detection"): search up to 3,000 nt upstream of the
RT start codon for short words recurring >=3 times at regular spacing (100-450 nt,
start to start). For the 95 ART members: the 20 most frequent 14-nt words of the
window serve as seeds (>=3 distinct bases, no homopolymer of six or more). Copies are
non-overlapping windows within two mismatches of a seed. An array is called when the
longest run of regularly spaced copies, R, is >=3 and exceeds the longest such run in
100 mononucleotide shuffles of the window. No array is called when a longer run of
copies spaced under 100 nt apart is present. R is the copy number. Loci with <1,500 nt
of contig upstream and no array are "not_assessed".

Step 2 (delimitation, "Array delimitation and measurement"): rescan up to 6,000 nt
upstream with every recurring 10-nt word as a seed; copies at <=1 mismatch; longest
chain at near-constant spacing (60-600 nt, 30% tolerance, single skipped copies
allowed); score = (copies - 1) x information content of the conserved block relative
to the region's base composition; retain when the score exceeds the best chain score
in each of 200 shuffles (50-nt block permutation, preserving local composition), in
both a 3,000-nt and a 6,000-nt window; chains scoring <2x the best shuffle are
retested against 2,000 shuffles.

INTERPRETED marks the few points where the paper's wording leaves an implementation
choice open; each carries its reading inline.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field

REPEAT_LEN_SCAN = 14  # 14-nt seed words
SEEDS_PER_WINDOW = 20
COPY_MISMATCHES = 2
MIN_SPACING, MAX_SPACING = 100, 450
MIN_RUN = 3
N_SHUFFLES_SCAN = 100
MIN_DISTINCT_BASES = 3
HOMOPOLYMER = 6
MIN_UPSTREAM_FOR_ASSESSMENT = 1_500
MAX_UPSTREAM_SCAN = 3_000
MAX_UPSTREAM_DELIMIT = 6_000
SEED_LEN_DELIMIT = 10
# paper Methods p.32: an array is adjacent to the RT when "no annotated
# gene of 300 nt or more lay between its last copy and the RT"
RT_ADJACENCY_GENE_MIN_NT = 300
DELIMIT_MISMATCHES = 1
DELIMIT_SPACING = (60, 600)
SPACING_TOLERANCE = 0.30
N_SHUFFLES_DELIMIT = 200
N_SHUFFLES_DELIMIT_RETEST = 2_000
CONSENSUS_FRACTION = 0.80

BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def hamming_within(a: str, b: str, max_mm: int) -> bool:
    """True when a and b differ in at most max_mm positions (early exit)."""
    mm = 0
    for x, y in zip(a, b, strict=False):
        if x != y:
            mm += 1
            if mm > max_mm:
                return False
    return True


def mononucleotide_shuffles(seq: str, n: int, rng: random.Random) -> list[str]:
    """Permute the window's nucleotides at random, preserving base composition."""
    out = []
    chars = list(seq)
    for _ in range(n):
        rng.shuffle(chars)
        out.append("".join(chars))
    return out


def block_shuffles(seq: str, n: int, block: int, rng: random.Random) -> list[str]:
    """Permute nucleotides within consecutive 50-nt blocks (local composition kept)."""
    out = []
    for _ in range(n):
        chars = list(seq)
        for start in range(0, len(chars), block):
            seg = chars[start:start + block]
            rng.shuffle(seg)
            chars[start:start + block] = seg
        out.append("".join(chars))
    return out


def _seed_ok(seed: str) -> bool:
    if len(set(seed)) < MIN_DISTINCT_BASES:
        return False
    run = 1
    for i in range(1, len(seed)):
        run = run + 1 if seed[i] == seed[i - 1] else 1
        if run >= HOMOPOLYMER:
            return False
    return True


def _most_frequent_words(window: str, word_len: int, top: int) -> list[str]:
    counts: dict[str, int] = {}
    for i in range(len(window) - word_len + 1):
        w = window[i:i + word_len]
        counts[w] = counts.get(w, 0) + 1
    words = [w for w, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
             if _seed_ok(w)]
    return words[:top]


def _copies_of(window: str, seed: str, max_mm: int) -> list[int]:
    """Greedy left-to-right non-overlapping copy positions."""
    pos: list[int] = []
    i = 0
    while i + len(seed) <= len(window):
        if hamming_within(window[i:i + len(seed)], seed, max_mm):
            pos.append(i)
            i += len(seed)  # non-overlapping
        else:
            i += 1
    return pos


def _tolerant_window(run: list[int]) -> list[int] | None:
    """Longest contiguous sub-window of ``run`` whose gaps are all within 30% of
    the window's median gap. A single irregular gap must not discard an otherwise
    regular run: 0,200,400,700,900,1100 holds two compliant 3-copy windows even
    though the whole run fails the median check (S1 finding B2)."""
    gaps = [b - a for a, b in zip(run, run[1:], strict=False)]
    med = _median(gaps)
    if all(abs(g - med) <= SPACING_TOLERANCE * med for g in gaps):
        return run
    for size in range(len(run) - 1, MIN_RUN - 1, -1):
        for s in range(len(run) - size + 1):
            sub = run[s:s + size]
            g = [b - a for a, b in zip(sub, sub[1:], strict=False)]
            m = _median(g)
            if all(abs(x - m) <= SPACING_TOLERANCE * m for x in g):
                return sub
    return None


def _median(xs: list[float]) -> float:
    """True median: with two gaps the upper-middle index rejected 3-copy runs
    whose gap ratio was in [0.54, 0.71] (S1 finding B4: R >= 3 is the paper's
    own call threshold, so legal arrays must not fail the tolerance rule)."""
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def _regular_runs(positions: list[int]) -> list[list[int]]:
    """Split sorted copy positions into maximal runs with spacing in
    [MIN_SPACING, MAX_SPACING]; each run is reduced to its longest window whose
    gaps are within 30% of its median gap (INTERPRETED: "regularly spaced" —
    tolerance taken from the paper's delimitation step, which fixes 30%)."""
    if len(positions) < MIN_RUN:
        return []
    runs: list[list[int]] = []
    cur = [positions[0]]
    for prev, nxt in zip(positions, positions[1:], strict=False):
        gap = nxt - prev
        if MIN_SPACING <= gap <= MAX_SPACING:
            cur.append(nxt)
        else:
            if len(cur) >= MIN_RUN:
                runs.append(cur)
            cur = [nxt]
    if len(cur) >= MIN_RUN:
        runs.append(cur)
    kept = []
    for run in runs:
        win = _tolerant_window(run)
        if win is not None and len(win) >= MIN_RUN:
            kept.append(win)
    return kept


@dataclass
class ScanCall:
    locus: str
    status: str  # "array" | "no_array" | "not_assessed"
    R: int = 0
    seed: str | None = None
    copies: list[int] = field(default_factory=list)
    note: str = ""


def kmer_scan(locus: str, upstream: str, rng: random.Random,
              word_len: int = REPEAT_LEN_SCAN) -> ScanCall:
    """Paper step 1. ``upstream`` is the window 5' of the RT start codon (already
    oriented so the RT is downstream); fewer than MIN_UPSTREAM_FOR_ASSESSMENT nt
    upstream with no array yields status "not_assessed". The scanned window is the
    3,000 nt immediately adjacent to the RT (grok review 2026-09-24: the prefix
    would scan the far end when upstream exceeds the cap)."""
    window = upstream[-MAX_UPSTREAM_SCAN:].upper() if len(upstream) > MAX_UPSTREAM_SCAN \
        else upstream.upper()
    best = ScanCall(locus=locus, status="no_array")
    for seed in _most_frequent_words(window, word_len, SEEDS_PER_WINDOW):
        positions = _copies_of(window, seed, COPY_MISMATCHES)
        runs = _regular_runs(positions)
        if not runs:
            continue
        long_run = max(runs, key=len)
        # suppression: a longer run of copies spaced under MIN_SPACING apart wins
        short_runs = _runs_with_spacing_at_most(positions, MIN_SPACING - 1)
        if short_runs and max(len(r) for r in short_runs) > len(long_run):
            continue
        r = len(long_run)
        if r < MIN_RUN:
            continue
        shuffled_best = max(
            (_longest_regular(_copies_of(s, seed, COPY_MISMATCHES))
             for s in mononucleotide_shuffles(window, N_SHUFFLES_SCAN, rng)),
            default=0,
        )
        if r >= MIN_RUN and r > shuffled_best and r > best.R:
            best = ScanCall(locus=locus, status="array", R=r, seed=seed,
                            copies=long_run)
    if best.status != "array" and len(upstream) < MIN_UPSTREAM_FOR_ASSESSMENT:
        best.status = "not_assessed"
    return best


def _runs_with_spacing_at_most(positions: list[int], max_gap: int) -> list[list[int]]:
    runs: list[list[int]] = []
    cur: list[int] = []
    for prev, nxt in zip(positions, positions[1:], strict=False):
        if nxt - prev <= max_gap:
            if not cur:
                cur = [prev]
            cur.append(nxt)
        else:
            if cur:
                runs.append(cur)
                cur = []
    if cur:
        runs.append(cur)
    return runs


def _longest_regular(positions: list[int]) -> int:
    runs = _regular_runs(positions)
    return max((len(r) for r in runs), default=0)


# --------------------------------------------------------------------------
# Step 2: delimitation
# --------------------------------------------------------------------------

@dataclass
class DelimitedArray:
    locus: str
    copy_starts: list[int]
    repeat: str  # consensus of the conserved block
    score: float
    shuffles_used: int
    spacings: list[int] = field(default_factory=list)
    block_offset: int = 0  # start of the conserved block within a copy
    # (grok round-2: PWM consumers must slice copy_start+block_offset, else they
    # align the seed prefix when s > 0)
    # paper Methods p.32: chains with all spacings multiples of 3 whose copies
    # lie within an annotated gene are "set aside as coding repeats"
    coding_repeat: bool = False
    # paper Methods p.32: adjacent to the RT when no annotated gene of
    # RT_ADJACENCY_GENE_MIN_NT nt or more lies between the last copy and the
    # RT; None = no gene annotation was supplied (rule not assessable)
    rt_adjacent: bool | None = None


def _information_content(columns: list[str], background: dict[str, float]) -> float:
    """Shannon IC of aligned columns against the region's base composition."""
    ic = 0.0
    log = math_log2
    n = len(columns)
    if n == 0:
        return 0.0
    for col in zip(*columns, strict=False):
        counts = {b: col.count(b) / n for b in BASES}
        ic += sum(p * log(p / background[b]) for b, p in counts.items() if p > 0)
    return ic


def math_log2(x: float) -> float:
    import math

    return math.log2(x)


def _consensus_block(copies_seqs: list[str]) -> tuple[int, int, str]:
    """Longest contiguous block of columns whose consensus base is carried by >=80%
    of copies, with one lapse tolerated. Returns (start, end_exclusive, consensus)."""
    if not copies_seqs:
        return 0, 0, ""
    width = min(len(s) for s in copies_seqs)
    cons = []
    for i in range(width):
        col = [s[i] for s in copies_seqs]
        # sorted() fixes tie order: max() keeps the first maximal
        # element, so modal-count ties resolve alphabetically, never by
        # set iteration order (PYTHONHASHSEED-dependent; audit finding 2)
        base = max(sorted(set(col)), key=col.count)
        cons.append((base, col.count(base) / len(col)))
    blocks: list[tuple[int, int]] = []
    start = None
    lapses = 0
    for i, (_, frac) in enumerate(cons):
        if frac >= CONSENSUS_FRACTION:
            if start is None:
                start = i
            # Paper p.32: "with one lapse tolerated" — ONE lapse per block
            # (a total budget), NOT per run of passing columns: a passing
            # column between two failures must not reset it.
        elif start is not None and lapses == 0:
            lapses += 1  # the one tolerated lapse
        else:
            if start is not None:
                blocks.append((start, i))  # end exclusive: the second failing
                # column is not part of the repeat (the ONE tolerated lapse is)
            start = None
            lapses = 0
    if start is not None:
        # NOT-IN-PAPER: whether a block whose LAST column is the tolerated
        # lapse (no second failure follows) includes that trailing lapse
        # column — the paper's wording does not settle it; included here.
        blocks.append((start, len(cons)))
    if not blocks:
        return 0, 0, ""
    s, e = max(blocks, key=lambda b: b[1] - b[0])
    return s, e, "".join(cons[i][0] for i in range(s, e))


def _extend_chain(positions: list[int], i: int, j: int, first_eff: float,
                  first_skipped: bool) -> list[int] | None:
    """Extend the seed pair (i, j) — whose first gap is already interpreted as
    ``first_eff`` (regular or halved) — into the longest compliant chain."""
    cur = [positions[i], positions[j]]
    gaps: list[float] = [first_eff]
    skipped = first_skipped
    for nxt in positions[j + 1:]:
        gap = nxt - cur[-1]
        med = _median(gaps)
        if abs(gap - med) <= SPACING_TOLERANCE * med:
            eff = gap  # regular spacing
        elif not skipped and abs(gap / 2 - med) <= SPACING_TOLERANCE * med:
            eff = gap / 2  # a single skipped copy: half-gap ~ median
            skipped = True
        else:
            break
        gaps.append(eff)
        cur.append(nxt)
    return cur if len(cur) >= MIN_RUN else None


def _chains(positions: list[int]) -> list[list[int]]:
    """Chains of copies at near-constant spacing (DELIMIT_SPACING, 30% tolerance
    against the running median gap, a single skipped copy allowed).

    The skipped-copy rule is anchored on the running median: with ART spacings of
    60-600 nt a missing copy produces a gap of ~2x the median, which is still
    inside [60, 600] — so a skip can never be detected from range membership alone
    (S1 finding B1: the old outside-range test was unreachable and the 30% check
    then broke every such chain).

    The FIRST gap has no median to anchor on, and its reading can be ambiguous:
    a 400-nt opening gap may be regular period-400 spacing or a period-200 chain
    with its first copy skipped (S3 finding 6: 0,400,600 must chain). Both
    readings are explored; the 30% check on the following gaps arbitrates.

    Chains consume CONSECUTIVE entries of ``positions`` (S4 finding 2: seeding
    a first gap from positions[j], j > i+1, jumped over DETECTED copies without
    counting them — a skipped copy is a missed detection, not a detected-but-
    ignored one)."""
    lo, hi = DELIMIT_SPACING
    out: dict[tuple[int, ...], list[int]] = {}
    for i in range(len(positions) - 1):
        j = i + 1  # consecutive detections only
        gap = positions[j] - positions[i]
        candidates: list[tuple[float, bool]] = []
        if lo <= gap <= hi:
            candidates.append((float(gap), False))
        if lo <= gap / 2 <= hi:
            candidates.append((gap / 2, True))  # provisional first-gap skip
        for eff, skipped in candidates:
            chain = _extend_chain(positions, i, j, eff, skipped)
            if chain is not None:
                out[tuple(chain)] = chain
    return list(out.values())


def _chain_score(window: str, chain: list[int], seed_len: int) -> tuple[float, str]:
    copies_seqs = [window[p:p + seed_len] for p in chain]
    background = {b: window.count(b) / max(1, len(window)) for b in BASES}
    s, e, _ = _consensus_block(copies_seqs)
    columns = [seq[s:e] for seq in copies_seqs]
    ic = _information_content(columns, background)
    return (len(chain) - 1) * ic, "".join(
        max(sorted(set(col)), key=col.count)
        for col in zip(*columns, strict=False))


def _delimitation_flags(
    chain: list[int],
    spacings: list[int],
    window_len: int,
    gene_spans: list[tuple[int, int]] | None,
    rt_offset: int | None,
) -> tuple[bool, bool | None]:
    """Both paper Methods p.32 delimitation rules for one chain.

    rule 1: "Chains whose spacings were all multiples of three and whose
    copies lay within an annotated gene were set aside as coding repeats."
    rule 2: "An array was called adjacent to the RT when no annotated gene of
    300 nt or more lay between its last copy and the RT."
    INTERPRETED (F3, devin review PR#18): "lay between" = the gene lies fully
    inside the open gap (start >= last copy end, end <= RT position); a gene
    containing the last copy or crossing the RT is not "between". Spans are
    half-open [a, b) in window coordinates.
    """
    coding_repeat = bool(gene_spans) and all(s % 3 == 0 for s in spacings) and any(
        a <= chain[0] and chain[-1] + SEED_LEN_DELIMIT <= b for a, b in gene_spans
    )
    rt_adjacent = None
    if gene_spans is not None:
        rt_pos = window_len if rt_offset is None else rt_offset
        last_copy_end = chain[-1] + SEED_LEN_DELIMIT
        rt_adjacent = not any(
            b - a >= RT_ADJACENCY_GENE_MIN_NT
            and a >= last_copy_end and b <= rt_pos
            for a, b in gene_spans
        )
    return coding_repeat, rt_adjacent


def delimit_array(
    locus: str,
    upstream: str,
    rng: random.Random,
    gene_spans: list[tuple[int, int]] | None = None,
    rt_offset: int | None = None,
) -> DelimitedArray | None:
    """Paper step 2. Tests every recurring 10-nt word as seed; retains the longest
    near-constant-spaced chain whose score beats the best chain in 200 (or 2,000 on
    weak margins) 50-nt-block shuffles, in both 3,000- and 6,000-nt windows.

    gene_spans: annotated gene (start, end) coordinates in the SAME frame as the
    returned copy_starts (i.e. positions in the upstream window, not genomic
    coordinates — the caller maps them). When given, both paper Methods p.32
    delimitation rules are applied: all-modulo-3 chains inside one annotated gene
    are flagged coding_repeat, and rt_adjacent reports whether any annotated gene
    of RT_ADJACENCY_GENE_MIN_NT nt or more lies between the last copy and the RT
    (rt_offset defaults to the end of the window — the RT immediately follows the
    supplied upstream). Without gene_spans both stay neutral (False / None)."""
    window6 = upstream[-MAX_UPSTREAM_DELIMIT:].upper() if len(upstream) > MAX_UPSTREAM_DELIMIT \
        else upstream.upper()
    window3 = window6[-MAX_UPSTREAM_SCAN:]
    counts: dict[str, int] = {}
    for i in range(len(window6) - SEED_LEN_DELIMIT + 1):
        w = window6[i:i + SEED_LEN_DELIMIT]
        counts[w] = counts.get(w, 0) + 1
    seeds = [w for w, c in counts.items() if c >= 2]  # INTERPRETED: "recurring" = >=2

    # Collect every candidate chain from every recurring seed; the paper retains
    # the LONGEST near-constant-spaced chain "for scoring" (Methods p.32) — ties
    # resolved by the higher information score (INTERPRETED: unspecified).
    candidates: list[list[int]] = []
    for seed in seeds:
        positions = _copies_of(window6, seed, DELIMIT_MISMATCHES)
        candidates.extend(_chains(positions))
    if not candidates:
        return None
    chain = max(candidates, key=lambda c: (len(c),
                                           _chain_score(window6, c, SEED_LEN_DELIMIT)[0]))
    score, _ = _chain_score(window6, chain, SEED_LEN_DELIMIT)
    if score <= 0:
        return None
    shuffles_used = N_SHUFFLES_DELIMIT

    def best_shuffle_score(window: str, n: int) -> float:
        """Null = each shuffled region's own best chain: seed discovery, chaining,
        and scoring are redone inside every shuffle (grok review 2026-09-24)."""
        top = 0.0
        for shuf in block_shuffles(window, n, 50, rng):
            w_counts: dict[str, int] = {}
            for i in range(len(shuf) - SEED_LEN_DELIMIT + 1):
                word = shuf[i:i + SEED_LEN_DELIMIT]
                w_counts[word] = w_counts.get(word, 0) + 1
            for seed in (w for w, c in w_counts.items() if c >= 2):
                positions = _copies_of(shuf, seed, DELIMIT_MISMATCHES)
                for cand in _chains(positions):
                    s_val, _ = _chain_score(shuf, cand, SEED_LEN_DELIMIT)
                    top = max(top, s_val)
        return top

    base_max3 = best_shuffle_score(window3, N_SHUFFLES_DELIMIT)
    base_max6 = best_shuffle_score(window6, N_SHUFFLES_DELIMIT)
    beat = score > base_max3 and score > base_max6
    null_margin = max(base_max3, base_max6)
    if beat and score < 2 * null_margin:
        beat = all(score > best_shuffle_score(w, N_SHUFFLES_DELIMIT_RETEST)
                   for w in (window3, window6))
        shuffles_used = N_SHUFFLES_DELIMIT_RETEST
    if not beat:
        return None

    copies_seqs = [window6[p:p + SEED_LEN_DELIMIT] for p in chain]
    blk_s, _, repeat = _consensus_block(copies_seqs)
    spacings = [b - a for a, b in zip(chain, chain[1:], strict=False)]
    coding_repeat, rt_adjacent = _delimitation_flags(
        chain, spacings, len(window6), gene_spans, rt_offset)
    return DelimitedArray(locus=locus, copy_starts=chain, repeat=repeat,
                          score=score, shuffles_used=shuffles_used,
                          spacings=spacings, block_offset=blk_s,
                          coding_repeat=coding_repeat, rt_adjacent=rt_adjacent)


# --------------------------------------------------------------------------
# PWM extension and cross-array grouping (paper p.32, end of delimitation)
# --------------------------------------------------------------------------

def build_pwm(alignments: list[str], background: dict[str, float],
              pseudocount: float = 1e-3) -> list[dict[str, float]]:
    """Log-odds PWM over aligned copy sequences (columns = positions).

    NOT-IN-PAPER: pseudocount value; the paper specifies only that a PWM of the
    repeat was built.
    """
    width = min(len(s) for s in alignments)
    pwm = []
    for i in range(width):
        col = {b: pseudocount for b in BASES}
        for s in alignments:
            col[s[i]] += 1.0
        total = sum(col.values())
        pwm.append({b: math_log2((col[b] / total)
                                 / max(background.get(b, 0.25), 1e-6))
                    for b in BASES})  # a 0-frequency background base must not
        # divide by zero (windows with no T are legal input)
    return pwm


def pwm_score(seq: str, pwm: list[dict[str, float]]) -> float:
    return sum(pwm[i][b] for i, b in enumerate(seq) if i < len(pwm))


def pwm_max_shuffle_score(window: str, pwm: list[dict[str, float]],
                          rng: random.Random, n: int = 200) -> float:
    """Threshold = best PWM match anywhere in each of 200 block-shuffled windows."""
    top = 0.0
    w = len(pwm)
    for shuf in block_shuffles(window, n, 50, rng):
        for i in range(len(shuf) - w + 1):
            top = max(top, pwm_score(shuf[i:i + w], pwm))
    return top


def pwm_extend(array: DelimitedArray, upstream: str, rng: random.Random,
               gene_spans: list[tuple[int, int]] | None = None,
               rt_offset: int | None = None,
               ) -> DelimitedArray:
    """'A position weight matrix of the repeat was then built, and matches that
    scored above the maximum of 200 shuffled regions were counted as copies'
    (Methods p.32). Returns a new DelimitedArray with extended copy list.

    Delimitation flags (coding_repeat / rt_adjacent) are RECOMPUTED on the
    extended copy list when gene_spans is supplied (extension can push the
    last copy past the old annotation boundary). Without annotation the old
    flags are stale the moment copies are added: rt_adjacent resets to None
    (unverifiable) rather than silently keeping a value computed for a chain
    that no longer exists (devin review PR#18, F1)."""
    window = upstream[-MAX_UPSTREAM_DELIMIT:].upper()
    background = {b: window.count(b) / max(1, len(window)) for b in BASES}
    w = len(array.repeat)
    off = array.block_offset
    aligned = [window[p + off:p + off + w] for p in array.copy_starts
               if p + off + w <= len(window)]
    if len(aligned) < MIN_RUN:
        return array
    pwm = build_pwm(aligned, background)
    threshold = pwm_max_shuffle_score(window, pwm, rng)
    # Non-overlapping scan, left to right. The PWM models the conserved block,
    # which sits at copy_start + block_offset; a hit at window position i is a
    # BLOCK start, so the copy start is i - off (grok round-2 finding on the
    # known-set: `known` holds block starts and must be tested directly).
    starts = sorted(array.copy_starts)
    known_block_starts = {p + off for p in starts}
    i = 0
    extra: list[int] = []
    while i + w <= len(window):
        if i in known_block_starts:
            i += w  # delimited copy: skip its block
            continue
        if pwm_score(window[i:i + w], pwm) > threshold:
            if i >= off:  # a block start before block_offset would place the
                # copy start before the window itself (negative coordinate)
                extra.append(i - off)  # record copy starts, not block starts
            i += w
        else:
            i += 1
    if not extra:
        return array
    merged = sorted(set(starts) | set(extra))
    spacings = [b - a for a, b in zip(merged, merged[1:], strict=False)]
    if gene_spans is not None:
        coding_repeat, rt_adjacent = _delimitation_flags(
            merged, spacings, len(window), gene_spans, rt_offset)
    elif extra:
        # ANY extension changes spacings: interstitial copies can break the
        # all-mod-3 / single-gene coverage, so coding_repeat is unverifiable
        # (devin re-review N2 — guarding only on a moved last copy kept a
        # stale True). rt_adjacent survives when the last copy did not move.
        coding_repeat = False
        rt_adjacent = None if merged[-1] > starts[-1] else array.rt_adjacent
    else:
        coding_repeat, rt_adjacent = array.coding_repeat, array.rt_adjacent
    return DelimitedArray(locus=array.locus, copy_starts=merged, repeat=array.repeat,
                          score=array.score, shuffles_used=array.shuffles_used,
                          spacings=spacings, block_offset=off,  # offset must survive
                          coding_repeat=coding_repeat, rt_adjacent=rt_adjacent)


def cross_scan(arrays: list[DelimitedArray], upstreams: dict[str, str],
               rng: random.Random) -> dict[str, set[str]]:
    """'Each locus was also scanned with the position weight matrix of every other
    array to group arrays that share a repeat' (Methods p.32).

    Returns locus -> set of other loci whose PWM matches that locus's upstream
    region above the shuffled threshold. Edges are symmetrized by union: A in
    groups[B] iff B in groups[A].
    """
    pwms: dict[str, list[dict[str, float]]] = {}
    for arr in arrays:
        window = upstreams.get(arr.locus, "")[-MAX_UPSTREAM_DELIMIT:].upper()
        w = len(arr.repeat)
        off = arr.block_offset
        aligned = [window[p + off:p + off + w] for p in arr.copy_starts
                   if p + off + w <= len(window)]
        if aligned:
            pwms[arr.locus] = build_pwm(aligned,
                                        {b: window.count(b) / max(1, len(window))
                                         for b in BASES})
    groups: dict[str, set[str]] = {a.locus: set() for a in arrays}
    for arr in arrays:
        window = upstreams.get(arr.locus, "")[-MAX_UPSTREAM_DELIMIT:].upper()
        if not window:
            continue
        for other_locus, pwm in pwms.items():
            if other_locus == arr.locus:
                continue
            w = len(pwm)
            if len(window) < w:
                continue
            threshold = pwm_max_shuffle_score(window, pwm, rng)
            best = max((pwm_score(window[i:i + w], pwm)
                        for i in range(len(window) - w + 1)), default=0.0)
            if best > threshold:
                groups[arr.locus].add(other_locus)
                groups[other_locus].add(arr.locus)  # symmetrize by union
    return groups


def exact_word_scan(locus: str, upstream: str, word_len: int = 12) -> ScanCall:
    """Paper's second scan setting (Methods p.32): an exact word (default 12 nt)
    recurring three times at regular spacing (100-450 nt, start to start) calls an
    array — no shuffle control, no mismatch allowance. Applies to phylogeny-tip
    loci and as the R=3 retention rule; callers reach it via
    scan_with_exact_word_fallback."""
    window = upstream[-MAX_UPSTREAM_SCAN:].upper()
    counts: dict[str, list[int]] = {}
    for i in range(len(window) - word_len + 1):
        counts.setdefault(window[i:i + word_len], []).append(i)
    best: ScanCall | None = None
    for word, positions in counts.items():
        runs = _regular_runs(positions)
        if not runs:
            continue
        run = max(runs, key=len)
        if len(run) >= MIN_RUN and (best is None or len(run) > best.R):
            best = ScanCall(locus=locus, status="array", R=len(run), seed=word,
                            copies=run,
                            note="exact-word rule, no shuffle control")
    return best or ScanCall(locus=locus, status="no_array")


def aligned_repeat_length(copies_seqs: list[str],
                          mafft_exe: str | None = None) -> int | None:
    """Second repeat-length measurement of paper Methods p.32: "The repeat
    length was measured both on ungapped copies and on a MAFFT alignment of
    the copies." Shells out to mafft --auto (an explicit mafft_exe, else
    $ARTHARNESS_MAFFT, else PATH); returns the alignment width after
    trimming all-gap TERMINAL columns, or None when no mafft is available
    (the caller keeps the ungapped measurement alone). NOT-IN-PAPER: the
    trimming convention; live mafft is exercised via recorded-fixture tests
    per repo convention, never in the offline suite.
    """
    import os
    import shutil
    import subprocess

    exe = mafft_exe or os.environ.get("ARTHARNESS_MAFFT") or shutil.which("mafft")
    if not exe:
        return None
    fasta = "".join(f">c{i}\n{s}\n" for i, s in enumerate(copies_seqs))
    try:
        proc = subprocess.run([str(exe), "--auto", "-"], input=fasta,
                              capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    seqs = ["".join(line.strip() for line in block.splitlines()[1:])
            for block in proc.stdout.split(">")[1:]]  # strip headers
    if not seqs:
        return None
    width = max(len(s) for s in seqs)
    seqs = [s.ljust(width, "-") for s in seqs]
    cols = range(width)
    keep = [i for i in cols
            if any(s[i] != "-" for s in seqs)]  # any-copy base keeps a column
    if not keep:
        return 0
    # trim only TERMINAL all-gap columns; internal gaps count toward the
    # aligned length (the aligned block can exceed the ungapped consensus)
    return keep[-1] - keep[0] + 1


def scan_with_exact_word_fallback(locus: str, upstream: str,
                                  rng: random.Random) -> ScanCall:
    """Default scan setting with the paper's exact-word retention (Methods
    p.32): "One locus with R = 3 that did not exceed its shuffles was
    retained because an exact 12-nt word recurred three times at such
    spacing." kmer_scan first; a locus it scores no_array is re-tested by
    exact_word_scan. The <1,500-nt not_assessed rule applies to the combined
    outcome ("loci with less than 1,500 nt of contig upstream of the RT and
    no array were recorded as not assessed")."""
    call = kmer_scan(locus, upstream, rng)
    if call.status == "array":
        return call
    fallback = exact_word_scan(locus, upstream)
    if fallback.status == "array":
        return fallback
    if len(upstream) < MIN_UPSTREAM_FOR_ASSESSMENT:
        return ScanCall(locus=locus, status="not_assessed")
    return fallback
