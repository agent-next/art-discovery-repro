import importlib.util
import json
from pathlib import Path


def _load():
    spec = importlib.util.spec_from_file_location(
        "forensics", Path(__file__).parent.parent / "experiments" / "forensics.py")
    mod = importlib.util.module_from_spec(spec)
    # register before exec: dataclasses resolves cls.__module__ through
    # sys.modules — an unregistered module makes @dataclass raise NoneType
    import sys
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_scan_records_finds_identifier(tmp_path: Path):
    mod = _load()
    (tmp_path / "records" / "t0001").mkdir(parents=True)
    (tmp_path / "records" / "t0001" / "summary.md").write_text(
        "examined contig L0050 from Logan; nothing unusual.")
    (tmp_path / "records" / "t0002").mkdir()
    (tmp_path / "records" / "t0002" / "summary.md").write_text("no relevant hits here.")
    ids = tmp_path / "ids.txt"
    ids.write_text("L0050\nMW218148.1\n")
    hits = mod.scan_records(tmp_path / "records", mod.load_identifiers(ids))
    assert hits == {"t0001": ["L0050"]}


def test_scan_transcripts_counts_dna_and_remarks(tmp_path: Path):
    mod = _load()
    dna = "ACGT" * 60  # 240 nt contiguous run (>= paper's 200-nt threshold)
    (tmp_path / "sessions").mkdir()
    tr = tmp_path / "sessions" / "s1.log"
    tr.write_text(f"retrieved flank of L0050: {dna}\n"
                  "I can see by eye a tandem repeat array!\n"
                  "weather in SF is nice\n")
    ids = tmp_path / "ids.txt"
    ids.write_text("L0050\n")
    out = mod.scan_transcripts(tmp_path / "sessions", mod.load_identifiers(ids))
    assert len(out) == 1
    assert out[0]["dna_runs_ge200nt"] == 1
    assert out[0]["repeat_remarks"] == 1
    assert "tandem repeat" in out[0]["repeat_remark_samples"][0]


def test_short_dna_not_counted(tmp_path: Path):
    mod = _load()
    tr = tmp_path / "s2.log"
    tr.write_text("L0050 flank: ACGT" * 10)  # only 40-nt runs
    ids = tmp_path / "ids.txt"
    ids.write_text("L0050\n")
    out = mod.scan_transcripts(tmp_path / ".", mod.load_identifiers(ids))
    assert out[0]["dna_runs_ge200nt"] == 0


def test_cli_end_to_end(tmp_path: Path):
    mod = _load()
    (tmp_path / "records" / "t0001").mkdir(parents=True)
    (tmp_path / "records" / "t0001" / "summary.md").write_text("hit MW248466.1 here")
    ids = tmp_path / "ids.txt"
    ids.write_text("MW248466.1\n")
    outp = tmp_path / "findings.json"
    import sys

    sys.argv = ["forensics.py", "--records", str(tmp_path / "records"),
                "--identifiers", str(ids), "--out", str(outp)]
    mod.main()
    data = json.loads(outp.read_text())
    assert data["record_hits"] == {"t0001": ["MW248466.1"]}
    assert data["identifiers_searched"] == 1


def test_same_line_remark_before_dna_not_after(tmp_path: Path):
    mod = _load()
    (tmp_path / "s").mkdir()
    dna = "ACGT" * 60
    # remark BEFORE the DNA on the same line -> not "after DNA"
    tr = tmp_path / "s" / "a.log"
    tr.write_text(f"tandem repeat suspected then {dna}\n")
    ids = tmp_path / "ids.txt"
    ids.write_text("L0050\n")
    assert mod.scan_transcripts(tmp_path / "s", mod.load_identifiers(ids)) == []


def test_scan_transcripts_remark_after_dna_counts(tmp_path: Path):
    mod = _load()
    (tmp_path / "s").mkdir()
    dna = "ACGT" * 60
    tr = tmp_path / "s" / "a.log"
    tr.write_text(f"flank of L0050: {dna}\ntandem repeat array!\n")
    ids = tmp_path / "ids.txt"
    ids.write_text("L0050\n")
    out = mod.scan_transcripts(tmp_path / "s", mod.load_identifiers(ids))
    assert out[0]["repeat_remarks_after_dna"] == 1


def test_scan_records_word_boundary(tmp_path: Path):
    mod = _load()
    (tmp_path / "records" / "t0001").mkdir(parents=True)
    (tmp_path / "records" / "t0001" / "summary.md").write_text("contig XL0050Y only")
    ids = tmp_path / "ids.txt"
    ids.write_text("L0050\n")
    hits = mod.scan_records(tmp_path / "records", mod.load_identifiers(ids))
    assert hits == {}  # bare-substring would falsely match


def test_remark_after_dna_on_earlier_line_with_own_dna(tmp_path: Path):
    # grok round-4 probe: remark's own line also carries DNA, but line 0 already
    # had a >=200-nt DNA run -> the remark must still count as after-DNA
    mod = _load()
    (tmp_path / "s").mkdir()
    dna = "ACGT" * 60
    tr = tmp_path / "s" / "a.log"
    tr.write_text(f"L0050 flank: {dna}\ntandem repeat suspected then {dna}\n")
    ids = tmp_path / "ids.txt"
    ids.write_text("L0050\n")
    out = mod.scan_transcripts(tmp_path / "s", mod.load_identifiers(ids))
    assert out[0]["dna_runs_ge200nt"] == 2
    assert out[0]["repeat_remarks_after_dna"] == 1


def test_typed_identifier_sets_rt_vs_contig(tmp_path):
    # Paper p.38: the forensics used 130 RT ids and 171 contig ids as
    # DISTINCT sets — "sessions naming a contig" are the event-walk subset.
    # A flat untyped set cannot express that (raw-fidelity audit MISMATCH).
    mod = _load()
    (tmp_path / "rt.txt").write_text("RTX_001\nRTX_002\n")
    (tmp_path / "contigs.txt").write_text("L0050\nMW218148.1\nON921432.1\n")
    sets = mod.load_identifier_sets(tmp_path / "rt.txt", tmp_path / "contigs.txt")
    assert len(sets.rt_ids) == 2 and len(sets.contig_ids) == 3
    assert sets.all == {"RTX_001", "RTX_002", "L0050", "MW218148.1",
                        "ON921432.1"}
    # paper numbers recorded where the rerun can compare against them
    assert sets.paper_counts == {"rt_ids": 130, "contig_ids": 171}
    # the event-walk gate: which named identifiers are contigs
    assert sets.contig_named({"L0050", "RTX_001"}) == {"L0050"}
    assert sets.contig_named({"RTX_001"}) == set()
