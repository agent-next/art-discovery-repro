"""scripts/scan_rt_loci.py: per-locus array scan on GenBank records (paper Methods p.32)."""

from __future__ import annotations

import importlib.util
import random
import sys
from pathlib import Path

from artharness.arrays import revcomp

_spec = importlib.util.spec_from_file_location(
    "scan_rt_loci", Path(__file__).parent.parent / "scripts" / "scan_rt_loci.py")
s = importlib.util.module_from_spec(_spec)
sys.modules["scan_rt_loci"] = s
_spec.loader.exec_module(s)


def _dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _world(strand: str) -> tuple[str, int, int]:
    """Plus-strand genome: 2.4 kb background, a 6-copy array (32-nt repeat, ~150-nt
    spacing), 300 nt of spacer, then a 900-nt RT gene. Returns (seq, rt_lo, rt_hi);
    for '-' the whole genome is reverse-complemented so the RT sits on the other strand."""
    rng = random.Random(7)
    repeat = _dna(rng, 32)
    arr = "".join(repeat + _dna(rng, 118 + 5 * i) for i in range(6))
    head = _dna(rng, 2400) + arr + _dna(rng, 300)
    rt = "ATG" + _dna(rng, 894) + "TAA"
    seq = head + rt + _dna(rng, 500)
    lo, hi = len(head) + 1, len(head) + len(rt)
    if strand == "-":
        n = len(seq)
        return revcomp(seq), n - hi + 1, n - lo + 1
    return seq, lo, hi


def _genbank(seq: str, lo: int, hi: int, strand: str, pid: str = "YP_TEST01.1",
             other: str = "1..300") -> str:
    loc = f"complement({lo}..{hi})" if strand == "-" else f"{lo}..{hi}"
    rows = [f"        1 {seq[i:i + 60]}" for i in range(0, len(seq), 60)]
    return (f"LOCUS       TEST                {len(seq)} bp    DNA     linear   PHG 01-JAN-2026\n"
            "VERSION     NC_TEST.1\n"
            "SOURCE      Testphage\n"
            "  ORGANISM  Testphage jumbo\n"
            "FEATURES             Location/Qualifiers\n"
            f"     CDS             {loc}\n"
            f'                     /protein_id="{pid}"\n'
            f"     CDS             {other}\n"
            '                     /protein_id="YP_OTHER.1"\n'
            "ORIGIN\n" + "\n".join(rows) + "\n//\n")


def test_parse_record_reads_cds_strand_and_sequence():
    seq, lo, hi = _world("-")
    rec = s.parse_record(_genbank(seq, lo, hi, "-"))
    assert rec.accession == "NC_TEST.1" and rec.organism == "Testphage jumbo"
    rt = next(c for c in rec.cds if c.protein_id == "YP_TEST01.1")
    assert (rt.lo, rt.hi, rt.strand) == (lo, hi, "-")
    assert rec.seq == seq


def test_upstream_window_is_oriented_so_the_rt_follows_it():
    for strand in "+-":
        seq, lo, hi = _world(strand)
        rec = s.parse_record(_genbank(seq, lo, hi, strand))
        rt = next(c for c in rec.cds if c.protein_id == "YP_TEST01.1")
        window, _ = s.upstream_window(rec, rt)
        plus_seq, plus_lo, _ = _world("+")
        assert window == plus_seq[max(0, plus_lo - 1 - 6000):plus_lo - 1]


def test_gene_spans_are_mapped_into_window_coordinates():
    seq, lo, hi = _world("+")
    rec = s.parse_record(_genbank(seq, lo, hi, "+"))
    rt = next(c for c in rec.cds if c.protein_id == "YP_TEST01.1")
    window, spans = s.upstream_window(rec, rt)
    assert spans == [(0, 300)] and len(window) == lo - 1
    seq_m, lo_m, hi_m = _world("-")
    other_lo = hi_m + 1000
    rec_m = s.parse_record(_genbank(seq_m, lo_m, hi_m, "-", other=f"{other_lo}..{other_lo + 299}"))
    rt_m = next(c for c in rec_m.cds if c.protein_id == "YP_TEST01.1")
    w_m, spans_m = s.upstream_window(rec_m, rt_m)
    # 1000 nt past the RT's 5' end on the reverse strand, read in the RT's orientation
    assert spans_m == [(len(w_m) - 1299, len(w_m) - 999)]


def test_planted_array_is_found_on_both_strands_and_a_bare_window_is_not(tmp_path: Path):
    for strand in "+-":
        seq, lo, hi = _world(strand)
        rec = s.parse_record(_genbank(seq, lo, hi, strand))
        rt = next(c for c in rec.cds if c.protein_id == "YP_TEST01.1")
        row = s.scan_locus(rec, rt)
        assert row["status"] == "array" and row["R"] >= 5, (strand, row)
        assert row["strand"] == strand
    bare = _dna(random.Random(3), 4000) + "ATG" + _dna(random.Random(4), 897)
    rec = s.parse_record(_genbank(bare, 4001, 4900, "+"))
    assert s.scan_locus(rec, next(c for c in rec.cds if c.protein_id == "YP_TEST01.1"))[
        "status"] == "no_array"


def test_cli_scans_only_requested_ids_from_a_gzip_flat_file(tmp_path: Path):
    import gzip

    seq, lo, hi = _world("+")
    gb = tmp_path / "g.gbff.gz"
    with gzip.open(gb, "wt") as fh:
        fh.write(_genbank(seq, lo, hi, "+", "YP_TEST01.1"))
        fh.write(_genbank(seq, lo, hi, "+", "YP_SKIPPED.1"))
    ids = tmp_path / "ids.txt"
    ids.write_text("YP_TEST01.1\nYP_MISSING.9\n")
    out = tmp_path / "out.tsv"
    assert s.main([str(ids), str(gb), "-o", str(out)]) == 0
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == s.COLUMNS and len(lines) == 2
    assert lines[1].startswith("YP_TEST01.1\tNC_TEST.1\tTestphage jumbo\t+")
