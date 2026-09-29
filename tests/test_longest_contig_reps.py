"""Family definition, cluster representatives (paper Methods p.31): "the member on
the longest contig was taken as the representative of each cluster".

Power: MMseqs' own representative (first column of its cluster table) differs from
the longest-contig member in every fixture below, so a pass-through of MMseqs'
choice fails these tests.
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "pipeline" / "art_family" / "longest_contig_reps.py"
FAMILY = REPO / "pipeline" / "art_family" / "family_definition.sh"

_spec = importlib.util.spec_from_file_location("longest_contig_reps", SCRIPT)
lcr = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = lcr
_spec.loader.exec_module(lcr)


CLUSTERS = {
    "mmseqs_rep_A": ["mmseqs_rep_A", "A2", "A3"],
    "mmseqs_rep_B": ["mmseqs_rep_B", "B2"],
    "single": ["single"],
}
LENGTHS = {
    "mmseqs_rep_A": 5_000, "A2": 90_000, "A3": 12_000,
    "mmseqs_rep_B": 300, "B2": 300,
    "single": 1,
}


def test_longest_contig_member_wins_over_mmseqs_representative():
    reps = lcr.pick_representatives(CLUSTERS, LENGTHS)
    assert reps == ["A2", "B2", "single"]
    assert all(r != k for r, k in zip(reps[:2], CLUSTERS, strict=False))


def test_contig_length_tie_goes_to_smallest_id():
    reps = lcr.pick_representatives({"z": ["z", "b", "m"]},
                                    {"z": 700, "b": 700, "m": 700})
    assert reps == ["b"]


def test_missing_contig_length_fails_loudly_and_names_members():
    with pytest.raises(KeyError, match=r"2 cluster members lack a contig length.*x1.*x2"):
        lcr.pick_representatives({"r": ["r", "x1", "x2"]}, {"r": 10})


def test_cli_end_to_end_on_mmseqs_style_tables(tmp_path: Path):
    clu = tmp_path / "clu_cluster.tsv"
    clu.write_text("".join(f"{rep}\t{m}\n" for rep, ms in CLUSTERS.items() for m in ms))
    lens = tmp_path / "contig_len.tsv"
    lens.write_text("# protein\tcontig_len\n"
                    + "".join(f"{p}\t{n}\n" for p, n in LENGTHS.items()))
    out = tmp_path / "reps.ids"
    proc = subprocess.run([sys.executable, str(SCRIPT), str(clu), str(lens), str(out)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert out.read_text().split() == ["A2", "B2", "single"]
    assert "3 clusters -> 3 representatives" in proc.stderr


def test_cli_usage_error_exit_code(tmp_path: Path):
    proc = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True)
    assert proc.returncode == 2


def _dry_run_lines() -> list[str]:
    proc = subprocess.run(["bash", str(FAMILY), "--dry-run"], capture_output=True,
                          text=True, cwd=REPO)
    assert proc.returncode == 0, proc.stderr
    return proc.stdout.splitlines()


def test_family_definition_feeds_longest_contig_reps_into_famsa():
    lines = _dry_run_lines()
    idx = {key: next(i for i, ln in enumerate(lines) if key in ln)
           for key in ("mmseqs easy-cluster", "longest_contig_reps.py",
                       "art_reps.ids", "famsa ")}
    assert (idx["mmseqs easy-cluster"] < idx["longest_contig_reps.py"]
            < idx["famsa "])
    famsa = lines[idx["famsa "]]
    assert "art_reps.faa" in famsa
    assert "_rep_seq" not in "\n".join(lines)  # MMseqs' own representatives not used


def test_seed_blastp_searches_genbank_not_the_metagenomic_fasta():
    blastp = next(ln for ln in _dry_run_lines() if ln.startswith("blastp "))
    assert "genbank_protein" in blastp
    assert "seqdb.faa" not in blastp  # a FASTA is not a BLAST database
