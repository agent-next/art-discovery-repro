"""Flag-fidelity oracles for the pipeline shell wrappers.

pipeline/*.sh are parameter-exact wrappers of the Anthropic ART paper's
Methods. These tests assert the required paper flags literally appear in each
script's text, and that ``bash <script> --dry-run`` exits 0 while echoing at
least one real tool command. Offline only: --dry-run never executes tools.
"""

import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CENSUS = REPO / "pipeline" / "census"
ART_FAMILY = REPO / "pipeline" / "art_family"
RNASEQ = REPO / "pipeline" / "rnaseq"
DB = REPO / "pipeline" / "db"

SHELL_SCRIPTS = [
    CENSUS / "01_hmmsearch.sh",
    CENSUS / "03_cluster.sh",
    ART_FAMILY / "family_definition.sh",
    ART_FAMILY / "phylogeny.sh",
    RNASEQ / "sa1_infection.sh",
    DB / "build_subset.sh",
]

# A tool token each script's --dry-run stdout must contain.
DRY_RUN_TOOL = {
    "01_hmmsearch.sh": "hmmsearch",
    "03_cluster.sh": "mmseqs easy-cluster",
    "family_definition.sh": "diamond blastp",
    "phylogeny.sh": "iqtree3",
    "sa1_infection.sh": "bowtie2",
    "build_subset.sh": "prodigal-gv",
}


def _text(path: Path) -> str:
    return path.read_text()


def _assert_flags(path: Path, required: list[str]) -> None:
    text = _text(path)
    missing = [flag for flag in required if flag not in text]
    assert not missing, f"{path.name} missing paper flags/params: {missing}"


def test_01_hmmsearch_flags():
    script = CENSUS / "01_hmmsearch.sh"
    _assert_flags(
        script,
        ["hmmsearch", "--noali", "-Z 8", "--domZ 8", "-E 0.01", "--domE 0.01",
         "--domtblout"],
    )
    text = _text(script)
    # paper: 52 profile HMMs; expected yield 3,391,714 candidate hits
    assert "52" in text
    assert "3,391,714" in text


def test_03_cluster_flags():
    script = CENSUS / "03_cluster.sh"
    _assert_flags(
        script,
        ["mmseqs easy-cluster", "--min-seq-id 0.5", "-c 0.8", "--cov-mode 0"],
    )
    # paper: expected yield 198,290 clusters
    assert "198,290" in _text(script)


def test_family_definition_flags():
    _assert_flags(
        ART_FAMILY / "family_definition.sh",
        [
            "QQM14740.1",        # MarsHill RT seed (GenBank)
            "mafft --auto",      # phage-RT profile alignment
            "hmmbuild",
            "1e-5",              # hmmsearch/diamond E-value ceiling
            "350",               # length filter min (aa)
            "900",               # length filter max (aa)
            "2,019",             # myRT reference RTs pooled
            "--very-sensitive",  # diamond blastp mode
            "140",               # edge bitscore floor
            "117",               # MarsHill connected component size
            "400",               # cluster-expansion min length (aa)
            "823",               # full-length RTs after expansion
            "--min-seq-id 0.9",  # mmseqs 90% clustering
            "230",               # resulting cluster count
            "famsa",             # representative alignment
            "FastTree -lg",      # LG-model tree
            "500",               # >=500 nt non-coding upstream
            "ART_01",
            "ART_95",
            "genomad",           # geNomad end-to-end
        ],
    )


def test_phylogeny_flags():
    script = ART_FAMILY / "phylogeny.sh"
    _assert_flags(
        script,
        [
            "774",          # sequence-set size (or see split below)
            "98",           # ART+NART
            "617",          # retrons
            "59",           # other RTs
            "RVT_1",        # Pfam RT-domain profile
            "+/- 40",       # envelope +/- 40 residues
            "--localpair",  # MAFFT L-INS-i
            "-gappyout",    # trimal
            "254",          # trimmed columns
            "-m MFP",       # ModelFinder (BIC)
            "-B 1000",      # 1,000 ultrafast bootstrap
            "-alrt 1000",   # 1,000 SH-aLRT replicates
            "midpoint",     # midpoint rooting
        ],
    )


def test_sa1_infection_flags():
    _assert_flags(
        RNASEQ / "sa1_infection.sh",
        [
            "PRJNA836150",               # BioProject
            "--length_required 30",      # fastp: paper states ONLY min length 30
            "--very-sensitive",          # bowtie2
            "-X 1000",                   # bowtie2 max insert
            "--no-unal",                 # bowtie2
            "MW218148.1",                # SA1 genome
            "NZ_CP059679.1",             # S. lentus chromosome
            "-q 10",                     # samtools MAPQ >= 10
            "-f 2",                      # samtools: properly paired
            "tlen <= 1500",              # samtools: template <= 1,500 nt
            "-p -s 2 -t CDS",            # featureCounts: fragments on read-2 strand
            "259",                       # features: 258 CDS + array RNA
        ],
    )


def test_build_subset_flags():
    _assert_flags(
        DB / "build_subset.sh",
        [
            "prodigal-gv",     # CDS calling, default params
            "8000",            # exclude >8,000-residue proteins
            "tantan",          # low-complexity masking
            ">=50%",           # >=50% low-complexity exclusion
            "--approx-id 90",  # cascade: 90% identity
            "--approx-id 70",  # cascade: 70% identity
            "--approx-id 50",  # cascade: 50% identity
            "--member-cover 80",  # >=80% coverage of the shorter member (percent)
            "--mutual-cover 80",  # 50% step: >=80% mutual coverage (percent)
            "80th",            # rep = closest to 80th pct of cluster length
        ],
    )


@pytest.mark.parametrize("script", SHELL_SCRIPTS, ids=lambda p: p.name)
def test_shell_script_dry_run(script: Path):
    """--dry-run exits 0 and echoes at least one real tool command."""
    proc = subprocess.run(
        ["bash", str(script), "--dry-run"],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=REPO,
    )
    assert proc.returncode == 0, f"{script.name} --dry-run failed: {proc.stderr}"
    assert proc.stdout.strip(), f"{script.name} --dry-run printed nothing"
    assert DRY_RUN_TOOL[script.name] in proc.stdout


def test_build_subset_annotated_input_retains_source_gene_calls(tmp_path):
    # Paper p.28: prodigal-gv only "in unannotated assemblies"; "the source
    # gene calls were retained for annotated assemblies". ANNOTATED_FAA must
    # bypass prodigal-gv in the dry-run plan; without it prodigal-gv runs.
    import os
    proc = subprocess.run(
        ["bash", str(DB / "build_subset.sh"), "--dry-run"],
        capture_output=True, text=True, check=True)
    assert "prodigal-gv" in proc.stdout
    env = dict(os.environ, ANNOTATED_FAA="inputs/annotated.faa")
    proc2 = subprocess.run(
        ["bash", str(DB / "build_subset.sh"), "--dry-run"],
        capture_output=True, text=True, check=True, env=env)
    assert "prodigal-gv" not in proc2.stdout
    assert "retain source gene calls" in proc2.stdout
    assert "inputs/annotated.faa" in proc2.stdout
