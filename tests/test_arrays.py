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


def _extendable_upstream() -> str:
    # 4 near-constant copies (60-nt step), a 200-nt gap (not near-constant ->
    # the delimit chain stops at 4), then 2 more repeat copies the PWM
    # extension must pick up beyond the gap.
    rng = random.Random(3)
    parts = []
    for _ in range(4):
        parts.append("ACGTACGTTA" + "".join(rng.choice("AT") for _ in range(50)))
    parts.append("".join(rng.choice("GC") for _ in range(200)))
    for _ in range(2):
        parts.append("ACGTACGTTA" + "".join(rng.choice("AT") for _ in range(50)))
    return "".join(parts)


def test_pwm_extend_recomputes_delimitation_flags_with_annotation():
    # Devin F1: pwm_extend used to rebuild DelimitedArray without
    # coding_repeat/rt_adjacent — the flags silently reset on the extension
    # path, guaranteed to defeat both p.32 rules once callers pass annotation.
    up = _extendable_upstream()
    # annotation: one gene covers ONLY the original 4 copies (0..189); the
    # extension pushes the last copy to >= 300 -> no single gene covers all.
    gene_orig = [(0, 190)]
    arr = _fast_delimit("L_f", up, gene_spans=gene_orig)
    assert arr is not None and len(arr.copy_starts) == 4  # chain stops at gap
    assert arr.coding_repeat  # precondition: all 4 copies inside the gene
    ext = arrays.pwm_extend(arr, up, random.Random(7),
                            gene_spans=gene_orig, rt_offset=len(up) + 50)
    assert len(ext.copy_starts) > len(arr.copy_starts)  # extension happened
    assert ext.coding_repeat is False  # copies now escape the annotated gene
    assert ext.rt_adjacent is not None


def test_pwm_extend_marks_flags_stale_when_annotation_absent():
    up = _extendable_upstream()
    arr = _fast_delimit("L_f", up)
    assert arr is not None and arr.rt_adjacent is None  # no annotation in
    ext = arrays.pwm_extend(arr, up, random.Random(7))
    assert len(ext.copy_starts) > len(arr.copy_starts)  # extension happened
    # copies were added beyond the old last copy -> old flags would be stale;
    # without annotation they are unverifiable, not silently kept
    assert ext.rt_adjacent is None


def test_pwm_extend_interstitial_extension_still_resets_coding_flag():
    # devin re-review N2: a copy inserted BETWEEN existing ones changes the
    # spacings (mod-3 / gene coverage may break) even though the last copy
    # never moves — coding_repeat must reset, rt_adjacent must survive.
    up = ("ACGTACGTTA" + "AT" * 5) * 10  # copies every 20 nt (unchainable:
    # spacing < 60, so only pwm_extend sees them all)
    arr = arrays.DelimitedArray(
        locus="L_n2", copy_starts=[0, 60, 120, 180], repeat="ACGTACGTTA",
        score=100.0, shuffles_used=200, spacings=[60, 60, 60], block_offset=0,
        coding_repeat=True, rt_adjacent=True)
    ext = arrays.pwm_extend(arr, up, random.Random(7))
    assert len(ext.copy_starts) > 4  # interstitial copies picked up
    assert ext.copy_starts[-1] == 180  # last copy unchanged
    assert ext.coding_repeat is False  # spacings changed -> unverifiable
    assert ext.rt_adjacent is True  # last copy unmoved -> still valid


def test_rt_adjacency_gene_containing_last_copy_does_not_block():
    # Devin F3: "lay between its last copy and the RT" — a >=300-nt gene that
    # CONTAINS the last copy is not BETWEEN the copy and the RT and must not
    # block adjacency (overlap-with-gap semantics wrongly blocked it).
    up = _planted_mod3_upstream()
    rt = len(up) + 500
    last_copy_end = len(up) - 50 + 10  # last copy block end (approx)
    spanning_gene = [(last_copy_end - 400, last_copy_end + 50)]  # 450 nt, contains copy
    arr = _fast_delimit("L_x", up, gene_spans=spanning_gene, rt_offset=rt)
    assert arr is not None and arr.rt_adjacent is True  # not BETWEEN


# ---------------------------------------------------------------------------
# Paper p.32 second scan setting: "an array was called either from annotated
# repeat copies ... or from an exact 12-nt word that recurred three times at
# such spacing" + the R=3 retention ("One locus with R = 3 that did not
# exceed its shuffles was retained because an exact 12-nt word recurred
# three times"). exact_word_scan existed dead; wire + test it (batch2).
# ---------------------------------------------------------------------------

def test_exact_word_scan_calls_array_on_three_exact_copies():
    word = "ACGTACGTACGT"
    up = (word + "A" * 88) * 3  # 100-nt start-to-start spacing
    call = arrays.exact_word_scan("L_w", up)
    assert call.status == "array" and call.R == 3
    assert call.seed == word


def test_exact_word_scan_rejects_two_copies_and_wrong_spacing():
    word = "ACGTACGTACGT"
    two = (word + "A" * 88) * 2 + "TTTTTTTTTTTT"
    assert arrays.exact_word_scan("L_w2", two).status == "no_array"
    # spacing far outside the regular-run bounds -> no array
    wide = word + "A" * 900 + word + "A" * 900 + word
    assert arrays.exact_word_scan("L_w3", wide).status == "no_array"


def test_scan_with_exact_word_fallback_retains_shuffle_failures():
    # the R=3 rule: a locus whose chain does not beat its shuffles is still
    # retained when an exact 12-nt word recurs three times at such spacing
    word = "ACGTACGTACGT"
    up = (word + "A" * 88) * 3
    call = arrays.scan_with_exact_word_fallback("L_f2", up, random.Random(0))
    assert call.status == "array" and call.R == 3
    # short upstream + nothing found stays not_assessed (paper: "loci with
    # less than 1,500 nt of contig upstream of the RT and no array were
    # recorded as not assessed")
    short = arrays.scan_with_exact_word_fallback(
        "L_s", "ACGT" * 10, random.Random(0))
    assert short.status == "not_assessed"


def test_aligned_repeat_length_via_mafft_stub(tmp_path):
    # Paper p.32: "The repeat length was measured both on ungapped copies and
    # on a MAFFT alignment of the copies." The MAFFT leg shells out; this
    # test exercises the plumbing with a stub binary emitting a recorded
    # alignment (repo convention: live tools tested via recorded fixtures,
    # never executed offline). Terminal all-gap columns are trimmed; internal
    # gap columns count (the aligned block can be longer than the ungapped
    # consensus). Returns None when no mafft is available.
    stub = tmp_path / "mafft"
    stub.write_text("#!/bin/sh\n"
                    "echo '>c1'\necho 'ACGTACGTT-A'\n"
                    "echo '>c2'\necho 'ACGTACGTTA-'\n"
                    "echo '>c3'\necho 'ACGTAC-TT-A'\n")
    stub.chmod(0o755)
    copies = ["ACGTACGTTA", "ACGTACGTTA", "ACGTACTTA"]
    n = arrays.aligned_repeat_length(copies, mafft_exe=stub)
    assert n == 11  # full stub alignment width (internal gaps count)
    # all-gap terminal columns are trimmed
    stub2 = tmp_path / "mafft2"
    stub2.write_text("#!/bin/sh\n"
                     "echo '>c1'\necho '--ACGT--'\n"
                     "echo '>c2'\necho '--ACGT--'\n")
    stub2.chmod(0.755 * 1000 if False else 0o755)
    assert arrays.aligned_repeat_length(["ACGT", "ACGT"], mafft_exe=stub2) == 4
    # no mafft anywhere -> None (caller keeps the ungapped measurement only)
    assert arrays.aligned_repeat_length(copies, mafft_exe=tmp_path / "nope") is None


def test_ungapped_and_aligned_measurements_pair():
    stub_len = arrays.aligned_repeat_length(["ACGTACGTTA"] * 3,
                                            mafft_exe=None)  # explicit None
    assert stub_len is None  # None exe == unavailable, not an error


def test_exact_word_retention_fires_only_when_kmer_scan_declines():
    # devin round-3: the earlier fixture's spacers let kmer_scan itself call
    # an array (14-mers spanning word+spacer recurred), so the fallback never
    # fired — vacuous. Random 4-base spacers per copy + a precondition
    # assertion make the fallback the deciding rule.
    word = "ACGTACGTACGT"  # 12 nt: below the 14-mer scan, at exact-word length
    rng = random.Random(7)
    ends = [("AA", "AC"), ("CG", "GT"), ("TT", "TA")]  # distinct 14-mer flanks
    parts = []
    for pre2, post2 in ends:
        parts.append(word)
        parts.append(pre2 + "".join(rng.choice("ACGT") for _ in range(84)) + post2)
    up = "".join(parts)
    pre = arrays.kmer_scan("L_pre", up, random.Random(0))
    assert pre.status != "array"  # precondition: the DEFAULT setting declines
    call = arrays.scan_with_exact_word_fallback("L_r3", up, random.Random(0))
    assert call.status == "array" and call.R == 3  # retention rule fires
