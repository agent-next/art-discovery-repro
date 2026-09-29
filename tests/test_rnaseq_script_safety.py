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
    lines = out.splitlines()
    rm = [ln for ln in lines if ln.startswith("rm -rf")]
    assert len(rm) == 2
    # plain FASTQ goes as soon as fastp has read it, before the aligner runs
    raw_at = lines.index(rm[0])
    assert "out/fastq/SRR19152328_1.fastq" in rm[0]
    assert any(ln.startswith("fastp") for ln in lines[:raw_at])
    assert not any("bowtie2 --very-sensitive" in ln for ln in lines[:raw_at])
    assert "out/bam/SRR19152328.bam" in rm[1] and "out/trim/SRR19152328_1.fq.gz" in rm[1]
    assert not any("sorted.bam" in r for r in rm)
    keep = subprocess.run(["bash", str(SCRIPT), "--dry-run", str(accs)], capture_output=True,
                          text=True, env={**env, "KEEP_INTERMEDIATES": "1"}).stdout
    assert "rm -rf" not in keep


def test_aligner_streams_into_the_filter_without_an_intermediate_sam(tmp_path: Path):
    accs = tmp_path / "acc.txt"
    accs.write_text("SRR19152328\n")
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run", str(accs)],
                         capture_output=True, text=True).stdout
    assert ".sam" not in out and "-S " not in out
    assert "bowtie2 --very-sensitive -X 1000 --no-unal" in out and "|" in out
    assert "samtools view -b -q 10 -f 2" in out


def test_ena_source_downloads_gz_and_skips_fasterq_dump(tmp_path: Path):
    accs = tmp_path / "acc.txt"
    accs.write_text("SRR19152328\n")
    env = {"PATH": "/usr/bin:/bin", "FASTQ_SOURCE": "ena", "OUTDIR": "out",
           "KEEP_INTERMEDIATES": "0"}
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run", str(accs)],
                         capture_output=True, text=True, env=env).stdout
    assert "fasterq-dump" not in out and "prefetch" not in out
    assert "filereport?accession=$1" in out and "bash SRR19152328 out/fastq" in out
    assert "-i out/fastq/SRR19152328_1.fastq.gz -I out/fastq/SRR19152328_2.fastq.gz" in out
    assert "out/fastq/SRR19152328_1.fastq.gz out/fastq/SRR19152328_2.fastq.gz" in [
        ln for ln in out.splitlines() if ln.startswith("rm -rf")][0]


def test_ena_source_resolves_urls_from_the_real_filereport_layout(tmp_path: Path):
    # ENA's filereport TSV with fields=fastq_ftp has TWO columns (run_accession, urls);
    # a first attempt cut the whole line and produced a malformed URL.
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    log = tmp_path / "curl.log"
    (bin_ / "curl").write_text(
        "#!/bin/bash\n"
        'case "$*" in\n'
        '  *filereport*) printf "run_accession\\tfastq_ftp\\nSRR19152328\\t'
        'ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/028/SRR19152328/SRR19152328_1.fastq.gz;'
        'ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/028/SRR19152328/SRR19152328_2.fastq.gz\\n" ;;\n'
        f'  *) echo "$*" >> {log}; f="";'
        ' while [ $# -gt 0 ]; do [ "$1" = -o ] && f="$2"; shift; done; : > "$f" ;;\n'
        "esac\n")
    (bin_ / "featureCounts").write_text(
        "#!/bin/bash\n"
        'while [ $# -gt 0 ]; do [ "$1" = -o ] && f="$2"; shift; done\n'
        'printf "Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\tx.bam\\ng1\\tc\\t1\\t9\\t+\\t9\\t5\\n" > "$f"\n')
    for tool in ("bowtie2-build", "fastp", "bowtie2", "samtools"):
        (bin_ / tool).write_text("#!/bin/bash\ncat > /dev/null 2>&1 < /dev/null; exit 0\n")
    for f in bin_.iterdir():
        f.chmod(0o755)
    for name in ("a.fna", "b.fna"):
        (tmp_path / name).write_text(">x\nACGT\n")
    accs = tmp_path / "acc.txt"
    accs.write_text("SRR19152328\n")
    (tmp_path / "out" / "fastq").mkdir(parents=True)
    (tmp_path / "out" / "trim").mkdir()
    (tmp_path / "out" / "bam").mkdir()
    env = {"PATH": f"{bin_}:/usr/bin:/bin", "FASTQ_SOURCE": "ena", "OUTDIR": str(tmp_path / "out"),
           "REF_SA1": str(tmp_path / "a.fna"), "REF_HOST": str(tmp_path / "b.fna"),
           "FEATURES": str(tmp_path / "f.gtf")}
    res = subprocess.run(["bash", str(SCRIPT), str(accs)], capture_output=True, text=True,
                         env=env)
    assert res.returncode == 0, res.stderr
    urls = [ln.split()[-1] for ln in log.read_text().splitlines()]
    assert urls == [
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/028/SRR19152328/SRR19152328_1.fastq.gz",
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/028/SRR19152328/SRR19152328_2.fastq.gz"]
