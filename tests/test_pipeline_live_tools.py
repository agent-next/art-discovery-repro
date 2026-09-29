"""Run the pipeline wrappers for real, on tiny synthetic inputs, against real tools.

The dry-run/flag tests prove the scripts print the paper's flags; these prove the flags
are accepted by the actual binaries and that the wrapped chain produces the expected
files (REVIEW.md bug class 2: invented flags). Each test skips when its tools are not on
PATH, so `make check` stays offline and tool-free; run them with the bio tools
installed (see TASKS.md).

Inputs are random proteins with planted families (fixed seed), so every result has an
independent oracle: the planted membership.
"""

import random
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
AA = "ACDEFGHIKLMNPQRSTVWY"


def need(*tools: str):
    missing = [t for t in tools if shutil.which(t) is None]
    return pytest.mark.skipif(bool(missing), reason=f"missing tools: {', '.join(missing)}")


def random_protein(rng: random.Random, n: int) -> str:
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def mutate(rng: random.Random, seq: str, rate: float) -> str:
    return "".join(rng.choice(AA) if i and rng.random() < rate else c
                   for i, c in enumerate(seq))


def write_fasta(path: Path, records: dict[str, str]) -> Path:
    path.write_text("".join(f">{k}\n{v}\n" for k, v in records.items()))
    return path


def planted(rng: random.Random, prefix: str, n: int, length: int, rate: float,
            ) -> dict[str, str]:
    root = random_protein(rng, length)
    return {f"{prefix}{i:02d}": mutate(rng, root, rate) for i in range(n)}


def write_family_hmm(tmp: Path, family: dict[str, str], name: str = "fam") -> Path:
    """Substitution-only members are aligned by construction: Stockholm needs no aligner."""
    sto = tmp / f"{name}.sto"
    sto.write_text("# STOCKHOLM 1.0\n" + "".join(f"{k} {v}\n" for k, v in family.items())
                   + "//\n")
    hmm = tmp / f"{name}.hmm"
    subprocess.run(["hmmbuild", "-n", name, str(hmm), str(sto)], check=True,
                   capture_output=True)
    return hmm


def sh(script: str, *args: str, env: dict | None = None, cwd: Path | None = None):
    import os
    return subprocess.run(["bash", str(REPO / "pipeline" / script), *args],
                          capture_output=True, text=True, cwd=cwd,
                          env={**os.environ, **(env or {})})


def domtbl_targets(path: Path) -> set[str]:
    return {ln.split()[0] for ln in path.read_text().splitlines()
            if ln.strip() and not ln.startswith("#")}


@need("hmmbuild", "hmmsearch")
def test_census_hmmsearch_recovers_the_planted_family_and_nothing_else(tmp_path: Path):
    rng = random.Random(1)
    fam = planted(rng, "member", 8, 250, 0.10)
    decoys = {f"decoy{i:02d}": random_protein(rng, 250) for i in range(20)}
    hmm = write_family_hmm(tmp_path, fam)
    db = write_fasta(tmp_path / "seqdb.faa", {**fam, **decoys})
    out = tmp_path / "hits.domtbl"
    proc = sh("census/01_hmmsearch.sh", str(hmm), str(db), str(out))
    assert proc.returncode == 0, proc.stderr
    assert domtbl_targets(out) == set(fam)


@need("mmseqs")
def test_census_mmseqs_cluster_separates_planted_families(tmp_path: Path):
    rng = random.Random(2)
    a = planted(rng, "a", 5, 300, 0.05)
    b = planted(rng, "b", 5, 300, 0.05)
    faa = write_fasta(tmp_path / "in.faa", {**a, **b})
    proc = sh("census/03_cluster.sh", str(faa), str(tmp_path / "clu"), str(tmp_path / "tmp"))
    assert proc.returncode == 0, proc.stderr
    clusters: dict[str, set[str]] = {}
    for ln in (tmp_path / "clu_cluster.tsv").read_text().splitlines():
        rep, member = ln.split("\t")
        clusters.setdefault(rep, set()).add(member)
    assert sorted(map(frozenset, clusters.values()), key=sorted) == sorted(
        [frozenset(a), frozenset(b)], key=sorted)


@need("hmmbuild", "hmmsearch", "seqkit", "mafft", "trimal", "iqtree3", "gotree")
def test_phylogeny_wrapper_produces_a_midpoint_rooted_tree_of_all_members(tmp_path: Path):
    rng = random.Random(3)
    fam = planted(rng, "rt", 8, 300, 0.15)
    hmm = write_family_hmm(tmp_path, fam, name="RVT_1")
    faa = write_fasta(tmp_path / "rt.faa", fam)
    out = tmp_path / "phylo"
    proc = sh("art_family/phylogeny.sh", str(faa), str(out), env={"RVT1_HMM": str(hmm)})
    assert proc.returncode == 0, proc.stderr
    tree = (out / "rt_tree.midpoint.nwk").read_text()
    assert all(name in tree for name in fam)
    assert tree.strip().endswith(";")


@need("diamond")
def test_diamond_cluster_flags_from_build_subset_are_real(tmp_path: Path):
    """The 90/70/50% cascade flags exactly as written in build_subset.sh (the class of
    bug where `--cov-mode` was invented for `diamond cluster`)."""
    import re
    text = (REPO / "pipeline" / "db" / "build_subset.sh").read_text()
    flag_sets = re.findall(r"run diamond cluster -d \S+ -o \S+ (--approx-id \d+ --\w+-cover \d+)",
                           text)
    assert len(flag_sets) == 3, flag_sets
    rng = random.Random(4)
    faa = write_fasta(tmp_path / "p.faa",
                      {**planted(rng, "x", 6, 250, 0.03), **planted(rng, "y", 6, 250, 0.03)})
    subprocess.run(["diamond", "makedb", "--in", str(faa), "-d", str(tmp_path / "p")],
                   check=True, capture_output=True)
    for i, flags in enumerate(flag_sets):
        proc = subprocess.run(
            ["diamond", "cluster", "-d", str(tmp_path / "p.dmnd"), "-o",
             str(tmp_path / f"c{i}.tsv"), *flags.split()],
            capture_output=True, text=True)
        assert proc.returncode == 0, f"{flags}: {proc.stderr}"
        assert (tmp_path / f"c{i}.tsv").stat().st_size > 0


@need("hmmbuild", "hmmsearch", "mafft", "diamond", "mmseqs", "seqkit", "famsa", "FastTree")
def test_family_definition_end_to_end_with_longest_contig_representatives(tmp_path: Path):
    """Whole wrapper, real tools. Only `blastp` (needs a GenBank BLAST db) and `genomad`
    (needs its database) are replaced by stubs that create their output files."""
    rng = random.Random(5)
    seed_id = "QQM14740.1"
    root = random_protein(rng, 420)
    members = {f"db{i:02d}": mutate(rng, root, 0.08) for i in range(12)}
    members[seed_id] = root
    decoys = {f"dec{i:02d}": random_protein(rng, 420) for i in range(20)}
    phage = {f"ph{i:02d}": mutate(rng, root, 0.20) for i in range(12)}
    myrt = {f"ref{i:02d}": random_protein(rng, 420) for i in range(6)}

    data = tmp_path / "data"
    data.mkdir()
    write_fasta(data / "seed.faa", {seed_id: root})
    write_fasta(data / "phage.faa", phage)
    write_fasta(data / "myrt.faa", myrt)
    write_fasta(data / "seqdb.faa", {**members, **decoys})
    (data / "genomes.fna").write_text(">g\nACGT\n")
    clust = data / "clust.tsv"
    clust.write_text("".join(f"{seed_id}\t{m}\n" for m in members))
    contig = data / "contig_len.tsv"
    # every clustered protein needs a contig length; phage RTs join the component too
    contig.write_text("".join(f"{m}\t{1000 + 37 * i}\n"
                              for i, m in enumerate(sorted({**members, **phage}))))
    stubs = tmp_path / "stubs"
    stubs.mkdir()
    for tool, out_flag in (("blastp", "-out"), ("genomad", None)):
        body = ('#!/usr/bin/env bash\nwhile [ $# -gt 0 ]; do '
                f'[ "$1" = "{out_flag}" ] && touch "$2"; shift; done\n'
                if out_flag else '#!/usr/bin/env bash\nmkdir -p "$3"\n')
        (stubs / tool).write_text(body)
        (stubs / tool).chmod(0o755)

    import os
    out = tmp_path / "out"
    env = {"PATH": f"{stubs}:{os.environ['PATH']}", "SEED_FAA": str(data / "seed.faa"),
           "SEED_ID": seed_id, "PHAGE_RT_FAA": str(data / "phage.faa"),
           "MYRT_REF_FAA": str(data / "myrt.faa"), "SEARCHDB_FAA": str(data / "seqdb.faa"),
           "CLUST_TSV": str(clust), "CONTIG_LEN_TSV": str(contig),
           "GENOMES_FNA": str(data / "genomes.fna"), "GENOMAD_DB": str(data / "gdb"),
           "GENBANK_BLASTDB": str(data / "gb"), "OUTDIR": str(out)}
    proc = sh("art_family/family_definition.sh", env=env, cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr[-2000:]

    component = {ln for ln in (out / "component117.ids").read_text().split()}
    assert component >= set(members)  # the planted family is one connected component
    assert not component & set(decoys) and not component & set(myrt)

    clusters: dict[str, list[str]] = {}
    for ln in (out / "art_823_clu_cluster.tsv").read_text().splitlines():
        rep, member = ln.split("\t")
        clusters.setdefault(rep, []).append(member)
    lengths = dict(ln.split("\t") for ln in contig.read_text().splitlines())
    expected = {min(ms, key=lambda m: (-int(lengths[m]), m)) for ms in clusters.values()}
    assert set((out / "art_reps.ids").read_text().split()) == expected
    tree = (out / "art_reps.nwk").read_text()
    assert all(m in tree for m in expected)
