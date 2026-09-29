import random

from artharness.arrays import DelimitedArray, build_pwm, cross_scan, pwm_extend

# Realistic ART-like repeat: paper says 15-49 nt with a ~15-nt palindromic core and
# less-conserved edges (Methods p.32/Results p.7). GAATTCCTTAAG is its own reverse
# complement; we flank it with per-copy random edges.
CORE = "GAATTCCTTAAG"
ALT_CORE = "TGCAAGTACGTT"  # same information content, unrelated sequence
REPEAT_LEN = 4 + len(CORE) + 4


def _repeat(rng: random.Random, core: str) -> str:
    return ("".join(rng.choice("ACGT") for _ in range(4)) + core
            + "".join(rng.choice("ACGT") for _ in range(4)))


def _window_with(rng: random.Random, core: str, n_copies: int = 14) -> tuple[str, list[int]]:
    copies = [20 + i * 136 for i in range(n_copies)]
    end = copies[-1] + REPEAT_LEN + 60
    chars = [rng.choice("ACGT") for _ in range(end)]
    for p in copies:
        word = _repeat(rng, core)
        for j, b in enumerate(word):
            if p + j < len(chars):
                chars[p + j] = b
    return "".join(chars), copies


def test_build_pwm_recovers_consensus_core():
    rng = random.Random(0)
    aligned = [_repeat(rng, CORE) for _ in range(10)]
    pwm = build_pwm(aligned, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
    consensus = "".join(max(col, key=col.get) for col in pwm)
    assert consensus[4:16] == CORE  # edges are random; core is conserved
    core_ic = sum(max(col.values()) for col in pwm[4:16])
    edge_ic = sum(max(col.values()) for col in pwm[:4])
    assert core_ic > edge_ic * 3  # core columns carry far more information


def test_pwm_extend_recovers_degenerate_copies():
    rng = random.Random(1)
    # plant extra copies inside spacer gaps (positions 40..156 etc. are free)
    extra = [80, 360]
    window, copies = _window_with(rng, CORE)
    chars = list(window)
    for p in extra:
        for j, b in enumerate(_repeat(rng, CORE)):
            if p + j < len(chars):
                chars[p + j] = b
    window = "".join(chars)
    arr = DelimitedArray(locus="L", copy_starts=copies, repeat="X" * REPEAT_LEN,
                         score=10.0, shuffles_used=200)
    out = pwm_extend(arr, window, rng)
    recovered = set(out.copy_starts) - set(copies)
    assert any(any(abs(e - x) <= 1 for e in recovered) for x in extra), (
        f"extra copies {extra} not recovered; got {sorted(recovered)}")


def test_pwm_extend_is_noop_without_extra_copies():
    rng = random.Random(2)
    window, copies = _window_with(rng, CORE)
    arr = DelimitedArray(locus="L", copy_starts=copies, repeat="X" * REPEAT_LEN,
                         score=10.0, shuffles_used=200)
    out = pwm_extend(arr, window, rng)
    assert set(out.copy_starts) == set(copies)


def test_cross_scan_groups_loci_sharing_a_repeat():
    win_a, copies_a = _window_with(random.Random(3), CORE)
    win_b, copies_b = _window_with(random.Random(4), CORE)
    win_c, copies_c = _window_with(random.Random(6), ALT_CORE)
    mk = lambda locus, starts: DelimitedArray(  # noqa: E731
        locus=locus, copy_starts=starts, repeat="X" * REPEAT_LEN, score=9.0,
        shuffles_used=200)
    arrays = [mk("A", copies_a), mk("B", copies_b), mk("C", copies_c)]
    groups = cross_scan(arrays, {"A": win_a, "B": win_b, "C": win_c},
                        random.Random(5))
    assert "B" in groups["A"] and "A" in groups["B"]  # same core -> grouped
    assert "C" not in groups["A"] and "C" not in groups["B"]
    assert "A" not in groups["C"]


def test_consensus_block_second_lapse_excluded():
    # grok round-2 probe: column fractions 1.0, 1.0, 0.5 (tolerated lapse),
    # 0.6 (second lapse) -> block must end exclusive of the 4th column
    from artharness.arrays import _consensus_block

    cols = [
        "AAAAAAAAAA",  # 1.0 A
        "AAAAAAAAAC",  # 0.9 A
        "AAAACCCCCC",  # 0.5 A  <- tolerated lapse
        "AAAACCCCCCC"[:11],  # 0.45.. use explicit below instead
    ]
    cols[3] = "AACCCCCCCC"[:10]
    s0, e0, cons = _consensus_block([c + "GG" for c in cols])
    assert (s0, e0) == (0, 3), (s0, e0)
    assert cons == "AAA"


def test_pwm_extend_uses_block_offset():
    # repeat conserved block starts at offset 4 within each copy: PWM alignment
    # must slice copy_start+offset, not the seed prefix (grok round-2 #4)
    rng = random.Random(9)
    copies = [20 + i * 136 for i in range(14)]
    end = copies[-1] + 24 + 60
    core = "GAATTCCTTAAG"
    chars = [rng.choice("ACGT") for _ in range(end)]
    for p in copies:
        word = "TTTT" + core + "GGGG"
        for j, b in enumerate(word):
            chars[p + j] = b
    window = "".join(chars)
    arr = DelimitedArray(locus="L", copy_starts=copies, repeat=core, score=10.0,
                         shuffles_used=200, block_offset=4)
    out = pwm_extend(arr, window, rng)
    assert set(out.copy_starts) == set(copies)  # nothing new, nothing lost


def test_cross_scan_with_block_offset():
    def build(seed):
        r = random.Random(seed)
        copies = [20 + i * 136 for i in range(14)]
        end = copies[-1] + 24 + 60
        chars = [r.choice("ACGT") for _ in range(end)]
        for p in copies:
            word = "TTTT" + core + "GGGG"
            for j, b in enumerate(word):
                chars[p + j] = b
        return "".join(chars), copies
    core = CORE
    wa, ca = build(3)
    wb, cb = build(4)
    arrays = [
        DelimitedArray(locus="A", copy_starts=ca, repeat=core, score=9.0,
                       shuffles_used=200, block_offset=4),
        DelimitedArray(locus="B", copy_starts=cb, repeat=core, score=9.0,
                       shuffles_used=200, block_offset=4),
    ]
    groups = cross_scan(arrays, {"A": wa, "B": wb}, random.Random(5))
    assert "B" in groups["A"] and "A" in groups["B"]


def test_pwm_extend_preserves_block_offset():
    # grok round-3 major: the extended DelimitedArray must keep block_offset
    rng = random.Random(11)
    copies = [20 + i * 136 for i in range(14)]
    end = copies[-1] + 24 + 60
    core = "GAATTCCTTAAG"
    chars = [rng.choice("ACGT") for _ in range(end)]
    for p in copies:
        for j, b in enumerate("TTTT" + core + "GGGG"):
            chars[p + j] = b
    # one extra conserved block in a spacer gap
    for j, b in enumerate("TTTT" + core + "GGGG"):
        chars[100 + j] = b
    window = "".join(chars)
    arr = DelimitedArray(locus="L", copy_starts=copies, repeat=core, score=10.0,
                         shuffles_used=200, block_offset=4)
    out = pwm_extend(arr, window, rng)
    assert len(out.copy_starts) > len(copies)
    assert out.block_offset == 4  # must survive the extension round-trip


def test_consensus_block_nonadjacent_second_lapse_ends_block():
    # Paper p.32 (verbatim): "at least 80% of copies carry the consensus base,
    # with ONE lapse tolerated" — a single lapse budget per block, not per run
    # of passing columns. A passing column between two failures must NOT reset
    # the budget; the block ends exclusive of the second failing column.
    # Column fractions (over the 4 copies): 1.0, 0.5 (lapse), 1.0 (pass),
    # 0.5 (second lapse), then 0.5 to end -> (0, 3) under a total budget.
    # (The reset-after-pass reading returned (0, 8): with an alternating
    # pattern every other failing column was re-tolerated. Verified by devin
    # review F5; the commit message's (0, 4) counterfactual was wrong.)
    from artharness.arrays import _consensus_block

    cols = [
        "AAAAAAAA",  # copy 1
        "ACACACAC",  # copy 2
        "AAAAAAAA",  # copy 3
        "ACACACAC",  # copy 4
    ]
    s0, e0, _cons = _consensus_block(cols)
    assert (s0, e0) == (0, 3), (s0, e0)
    # (the lapse column's consensus letter is a tied count; determinism of
    # THAT letter is pinned separately by the tie-break test below)
