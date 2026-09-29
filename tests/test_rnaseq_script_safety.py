"""S3b regression: pipeline/rnaseq/sa1_infection.sh must not eval, and its
accession gate must reject command-injection payloads before any run() call.

Power demo (REVIEW.md bars): the pre-fix version of this script DID execute
`$(...)` payloads planted in the accessions file; these tests fail on that
version (eval present / gate absent) and pass on the fixed one.
"""

import re
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).parent.parent / "pipeline" / "rnaseq" / "sa1_infection.sh"


def test_no_eval_remains_in_script():
    text = SCRIPT.read_text()
    # the word may appear only inside the SECURITY comment, never as a command
    code_lines = [ln for ln in text.splitlines()
                        if not ln.strip().startswith("#")]
    assert not any("eval" in ln for ln in code_lines), "eval must be gone from code"


def test_dry_run_smoke_still_prints_paper_commands():
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run", "/dev/null"],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert "bowtie2 --very-sensitive -X 1000 --no-unal" in out.stdout
    # paper counting model: properly paired + template cap + fragment/read-2 strand
    assert "samtools view -b -q 10 -f 2" in out.stdout
    assert "tlen <= 1500" in out.stdout
    assert "featureCounts -p -s 2 -t CDS" in out.stdout
    assert "samtools sort -o" in out.stdout
    # the small-RNA fastp settings must NOT leak into the infection runs
    assert "--disable_adapter_trimming" not in out.stdout
    assert "--trim_poly_g" not in out.stdout


def test_accession_gate_rejects_injection_payload():
    # extract acc_ok() from the script and fire a real payload at it
    text = SCRIPT.read_text()
    m = re.search(r"acc_ok\(\) \{.*?\n\}", text, re.DOTALL)
    assert m, "acc_ok gate must exist"
    payload = "SRR$(touch PWNED_ART_TEST)"
    probe = Path("/tmp") / "PWNED_ART_TEST"
    if probe.exists():
        probe.unlink()
    bash = f"{m.group(0)}\nif acc_ok '{payload}'; then echo ACCEPTED; else echo REJECTED; fi"
    out = subprocess.run(["bash", "-c", bash], capture_output=True, text=True)
    assert "REJECTED" in out.stdout
    assert not probe.exists(), "payload executed -- gate failed"


def test_accession_gate_accepts_real_accessions():
    text = SCRIPT.read_text()
    m = re.search(r"acc_ok\(\) \{.*?\n\}", text, re.DOTALL)
    bash = (f"{m.group(0)}\n"
            "for acc in SRR1234567 ERR987654; do acc_ok \"$acc\" || exit 1; done; echo OK")
    out = subprocess.run(["bash", "-c", bash], capture_output=True, text=True)
    assert "OK" in out.stdout


def test_keep_intermediates_zero_deletes_only_this_librarys_scratch(tmp_path: Path):
    accs = tmp_path / "acc.txt"
    accs.write_text("SRR19152328\n")
    env = {"PATH": "/usr/bin:/bin", "KEEP_INTERMEDIATES": "0", "OUTDIR": "out"}
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run", str(accs)],
                         capture_output=True, text=True, env=env).stdout
    rm = [ln for ln in out.splitlines() if ln.startswith("rm -rf")]
    assert len(rm) == 1
    assert "out/bam/SRR19152328.sam" in rm[0] and "out/trim/SRR19152328_1.fq.gz" in rm[0]
    assert "sorted.bam" not in rm[0]
    keep = subprocess.run(["bash", str(SCRIPT), "--dry-run", str(accs)], capture_output=True,
                          text=True, env={**env, "KEEP_INTERMEDIATES": "1"}).stdout
    assert "rm -rf" not in keep
