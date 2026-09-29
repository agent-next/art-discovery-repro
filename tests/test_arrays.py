import random

import pytest

from artharness import arrays
from artharness.arrays import delimit_array, kmer_scan, mononucleotide_shuffles

REPEAT = "CATGTGTATCGCATGT"  # 16-nt repeat observed by the worker in task t0062
SPACER = 120  # within the paper's 100-450 nt regular-spacing window


def synth_array_window(copies: int = 14, mismatches: int = 1, spacer: int = SPACER,
                       seed: int = 7) -> str:
    rng = random.Random(seed)
    parts = []
    for _i in range(copies):
        seq = list(REPEAT)
        for pos in rng.sample(range(len(REPEAT)), mismatches):
            seq[pos] = rng.choice([b for b in "ACGT" if b != seq[pos]])
        parts.append("".join(seq))
        parts.append("".join(rng.choice("ACGT") for _ in range(spacer)))
    return "".join(parts)


def test_kmer_scan_calls_array_on_synthetic_locus():
    window = synth_array_window()
    call = kmer_scan("L_test", window, random.Random(0))
    assert call.status == "array"
    assert call.R == 14
    assert call.seed is not None and len(call.seed) == 14
    assert len(call.copies) == 14


def test_kmer_scan_rejects_shuffled_window():
    window = synth_array_window()
    shuffled = mononucleotide_shuffles(window, 1, random.Random(1))[0]
    call = kmer_scan("L_shuf", shuffled, random.Random(0))
    assert call.status == "no_array"


def test_short_upstream_is_not_assessed():
    window = synth_array_window(copies=2)  # too few copies to call
    call = kmer_scan("L_short", window[:900], random.Random(0))
    assert call.status == "not_assessed"


def test_homopolymer_and_low_complexity_seeds_excluded():
    from artharness.arrays import _seed_ok

    assert not _seed_ok("AAAAAAAAAAAAAA")
    assert not _seed_ok("ACACACACACACAC")  # only 2 distinct bases
    assert _seed_ok(REPEAT[:14])


def test_delimit_recovers_copies_and_repeat():
    rng = random.Random(2)
    window = synth_array_window(mismatches=1)
    # shrink the shuffle counts for test speed (constants are module-level on purpose)
    orig = (arrays.N_SHUFFLES_DELIMIT, arrays.N_SHUFFLES_DELIMIT_RETEST)
    arrays.N_SHUFFLES_DELIMIT, arrays.N_SHUFFLES_DELIMIT_RETEST = 4, 8
    try:
        arr = delimit_array("L_test", window, rng)
    finally:
        arrays.N_SHUFFLES_DELIMIT, arrays.N_SHUFFLES_DELIMIT_RETEST = orig
    assert arr is not None
    assert len(arr.copy_starts) == 14
    assert all(60 <= s <= 600 for s in arr.spacings)
    # consensus repeat must be a 10-mer close to some 10-mer of the planted repeat
    # (any of the 7 recurring seed offsets, not necessarily offset 0)
    assert len(arr.repeat) == 10
    best_mm = min(sum(a != b for a, b in zip(arr.repeat, REPEAT[o:o + 10], strict=False))
                  for o in range(7))
    assert best_mm <= 1


def test_delimit_rejects_random_window():
    rng = random.Random(3)
    window = "".join(rng.choice("ACGT") for _ in range(6000))
    orig = (arrays.N_SHUFFLES_DELIMIT, arrays.N_SHUFFLES_DELIMIT_RETEST)
    arrays.N_SHUFFLES_DELIMIT, arrays.N_SHUFFLES_DELIMIT_RETEST = 4, 8
    try:
        assert delimit_array("L_rand", window, rng) is None
    finally:
        arrays.N_SHUFFLES_DELIMIT, arrays.N_SHUFFLES_DELIMIT_RETEST = orig


def test_paper_anchor_sequence_is_found():
    # the t0062 worker read a 2,900-nt flank of locus L0050 whose copies carry
    # <=1 mismatch against CATGTGT[AT]TCGCATGT (paper p.30-31, Methods)
    flank = synth_array_window(copies=14, mismatches=1, spacer=117, seed=11)[:2900]
    call = kmer_scan("L0050", flank, random.Random(4))
    assert call.status == "array"
    assert pytest.approx(call.R, abs=2) == 14


# ---------------------------------------------------------------------------
# Paper Methods p.32 delimitation rules (raw-fidelity audit 2026-09-29):
# "Chains whose spacings were all multiples of three and whose copies lay
# within an annotated gene were set aside as coding repeats. An array was
# called adjacent to the RT when no annotated gene of 300 nt or more lay
# between its last copy and the RT."
# ---------------------------------------------------------------------------

def _planted_mod3_upstream(copies: int = 5, spacing: int = 60) -> str:
    rng = random.Random(7)
    parts = []
    for _ in range(copies):
        parts.append("ACGTACGTTA")
        parts.append("".join(rng.choice("AT") for _ in range(spacing - 10)))
    return "".join(parts)


def _fast_delimit(locus: str, upstream: str, **kwargs):
    orig = (arrays.N_SHUFFLES_DELIMIT, arrays.N_SHUFFLES_DELIMIT_RETEST)
    arrays.N_SHUFFLES_DELIMIT, arrays.N_SHUFFLES_DELIMIT_RETEST = 4, 8
    try:
        return delimit_array(locus, upstream, random.Random(7), **kwargs)
    finally:
        arrays.N_SHUFFLES_DELIMIT, arrays.N_SHUFFLES_DELIMIT_RETEST = orig


def test_coding_repeat_rule_flags_all_mod3_chain_inside_gene():
    up = _planted_mod3_upstream()  # spacing 60 -> all multiples of 3
    arr = _fast_delimit("L_x", up)
    assert arr is not None
    # no annotation -> flags stay neutral
    assert arr.coding_repeat is False
    assert arr.rt_adjacent is None
    # gene covering every copy -> set aside as coding repeat
    covering = [(0, len(up))]
    arr2 = _fast_delimit("L_x", up, gene_spans=covering)
    assert arr2 is not None and arr2.coding_repeat is True
    # gene NOT covering the copies -> not a coding repeat
    off_target = [(len(up) - 5, len(up) + 50)]
    arr3 = _fast_delimit("L_x", up, gene_spans=off_target)
    assert arr3 is not None and arr3.coding_repeat is False


def test_rt_adjacency_rule_requires_no_300nt_gene_in_gap():
    up = _planted_mod3_upstream()
    rt = len(up) + 500  # RT 500 nt downstream of the window
    # a >=300-nt annotated gene in the gap -> NOT adjacent
    big_gene = [(len(up) + 50, len(up) + 400)]  # 350 nt
    arr = _fast_delimit("L_x", up, gene_spans=big_gene, rt_offset=rt)
    assert arr is not None and arr.rt_adjacent is False
    # only short genes in the gap -> adjacent
    small_gene = [(len(up) + 50, len(up) + 200)]  # 150 nt
    arr2 = _fast_delimit("L_x", up, gene_spans=small_gene, rt_offset=rt)
    assert arr2 is not None and arr2.rt_adjacent is True
    # empty annotation but gene list given -> no gene between -> adjacent
    arr3 = _fast_delimit("L_x", up, gene_spans=[], rt_offset=rt)
    assert arr3 is not None and arr3.rt_adjacent is True
