"""pipeline/lib/run.sh: values reach commands as arguments, never as shell text.

Power demo (REVIEW.md bars): the old `run()` did `eval "$1"` on an interpolated
string, so each payload below executed; these tests fail on that helper.
"""

import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LIB = REPO / "pipeline" / "lib" / "run.sh"
EVAL_FREE = [
    REPO / "pipeline" / "art_family" / "family_definition.sh",
    REPO / "pipeline" / "art_family" / "phylogeny.sh",
    REPO / "pipeline" / "db" / "build_subset.sh",
    REPO / "pipeline" / "rnaseq" / "sa1_infection.sh",
    LIB,
]


def _bash(snippet: str, cwd: Path, dry_run: int = 0) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", "-c", f"set -euo pipefail; DRY_RUN={dry_run}; source {LIB}\n{snippet}"],
        cwd=cwd, capture_output=True, text=True)


def test_no_eval_in_pipeline_code():
    for script in EVAL_FREE:
        code = [ln for ln in script.read_text().splitlines()
                if not ln.strip().startswith("#")]
        assert not any(" eval " in f" {ln} " or ln.strip().startswith("eval ")
                       for ln in code), f"eval found in {script.name}"


def test_run_treats_payload_argument_as_data(tmp_path: Path):
    payload = "x; touch PWNED_RUN; $(touch PWNED_SUB)"
    out = _bash(f"run touch '{payload}'", tmp_path)
    assert out.returncode == 0, out.stderr
    assert (tmp_path / payload).exists()  # created as ONE file with that literal name
    assert not (tmp_path / "PWNED_RUN").exists()
    assert not (tmp_path / "PWNED_SUB").exists()


def test_run_sh_passes_values_only_as_positional_args(tmp_path: Path):
    src = tmp_path / "in.txt"
    src.write_text("hello\n")
    payload_out = "o; touch PWNED_SH"
    out = _bash(f"run_sh 'cat \"$1\" > \"$2\"' '{src}' '{payload_out}'", tmp_path)
    assert out.returncode == 0, out.stderr
    assert (tmp_path / payload_out).read_text() == "hello\n"
    assert not (tmp_path / "PWNED_SH").exists()


def test_run_sh_dry_run_shows_substituted_command_without_executing(tmp_path: Path):
    out = _bash('run_sh \'cat "$1" > "$2"\' in.faa out.faa', tmp_path, dry_run=1)
    assert out.stdout.strip() == 'cat "in.faa" > "out.faa"'
    assert not (tmp_path / "out.faa").exists()


def test_run_sh_pipeline_and_quoted_awk_run_for_real(tmp_path: Path):
    (tmp_path / "t.tbl").write_text("#c\nA x\nB y\nA z\n")
    script = "grep -v '^#' \"$1\" | awk '{print $1}' | sort -u > \"$2\""
    out = _bash(f"run_sh \"$(cat <<'EOS'\n{script}\nEOS\n)\" t.tbl ids.txt", tmp_path)
    assert out.returncode == 0, out.stderr
    assert (tmp_path / "ids.txt").read_text() == "A\nB\n"
