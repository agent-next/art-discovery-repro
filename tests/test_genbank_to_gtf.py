"""genbank_to_gtf.py builds the RNA-seq feature file (paper p.36-37: 258 CDS + array RNA)."""

import importlib.util
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location(
    "gb2gtf", Path(__file__).parent.parent / "pipeline" / "rnaseq" / "genbank_to_gtf.py")
g = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(g)

# Verbatim feature-table layout of GenBank MW218148.1 (SA1); sequences omitted.
GB = """LOCUS       MW218148              260727 bp    DNA     linear   PHG 03-FEB-2021
VERSION     MW218148.1
FEATURES             Location/Qualifiers
     source          1..260727
                     /organism="Staphylococcus phage vB_StaM_SA1"
     CDS             90..5312
                     /codon_start=1
                     /product="hypothetical protein"
                     /protein_id="QPI16920.1"
                     /translation="MVLMYLNEEETVYIEDYLYKPIEGKHTLVIKNDRNTPIEVNINR
                     KVESLKINKVYNHKVNNNSSYEFDFPIDDTIPPDILSLFIDKQNTTYTNSNDVLNVPY"
     CDS             complement(64666..64974)
                     /codon_start=1
                     /protein_id="QPI16994.1"
     CDS             complement(<86973..87404)
                     /locus_tag="SA1_0099"
ORIGIN
"""


def test_cds_coordinates_strand_and_ids():
    assert g.parse_cds(GB) == [("QPI16920.1", 90, 5312, "+"),
                               ("QPI16994.1", 64666, 64974, "-"),
                               ("SA1_0099", 86973, 87404, "-")]


def test_gtf_has_one_typed_row_per_cds_plus_each_array_strand():
    gtf = g.to_gtf(GB, [(9597, 10371, "+"), (9597, 10371, "-")]).splitlines()
    assert len(gtf) == 5
    assert gtf[0].split("\t")[:8] == ["MW218148.1", "art-harness", "CDS", "90", "5312", ".",
                                      "+", "0"]
    assert 'gene_id "ART_array_RNA_fwd"' in gtf[3] and 'gene_id "ART_array_RNA_rev"' in gtf[4]


def test_joined_locations_are_refused_rather_than_guessed():
    with pytest.raises(ValueError, match="unsupported CDS location"):
        g.parse_cds("     CDS             join(1..10,20..30)\n")


def test_cli_writes_the_file(tmp_path: Path):
    gb, out = tmp_path / "x.gb", tmp_path / "x.gtf"
    gb.write_text(GB)
    assert g.main([str(gb), str(out), "--array", "1", "5", "+"]) == 0
    assert out.read_text().count("\n") == 4
