"""A toy genome world with a known answer, so the demo can be scored.

The paper's world is 1.9 billion protein clusters; this one is nine synthetic
"phage" contigs. Planted in them are the same kinds of objects the paper hunts:

- a reverse-transcriptase (RT) gene, marked by the catalytic motif ``[YF][ADV]DD``
  (NOT-IN-PAPER: the paper uses 52 profile HMMs; a motif is the smallest stand-in
  that keeps the task shape "find the RT genes");
- upstream of some RTs, a tandem repeat array with regular spacing (the ART
  signature: paper Results p.4, Methods p.32);
- decoys that a careless analysis would call: a too-short RT fragment, an RT with
  too little upstream sequence to assess, and a dense tandem repeat whose spacing is
  under 100 nt (the paper's suppression rule must reject it).

The truth table lives outside the agents' sandbox; only the scorecard reads it.
"""

from __future__ import annotations

import json
import random
from dataclasses import asdict, dataclass
from pathlib import Path

_CODON = {
    "A": ["GCT", "GCC", "GCA"], "R": ["CGT", "CGC", "AGA"], "N": ["AAT", "AAC"],
    "D": ["GAT", "GAC"], "C": ["TGT", "TGC"], "Q": ["CAA", "CAG"], "E": ["GAA", "GAG"],
    "G": ["GGT", "GGC", "GGA"], "H": ["CAT", "CAC"], "I": ["ATT", "ATC", "ATA"],
    "L": ["CTT", "CTC", "TTA", "TTG"], "K": ["AAA", "AAG"], "M": ["ATG"],
    "F": ["TTT", "TTC"], "P": ["CCT", "CCC", "CCA"], "S": ["TCT", "TCC", "AGT"],
    "T": ["ACT", "ACC", "ACA"], "W": ["TGG"], "Y": ["TAT", "TAC"],
    "V": ["GTT", "GTC", "GTA"],
}
_AA = "ARNDCQEGHILKMFPSTWYV"
_AA_NO_MOTIF_CHARS = [a for a in _AA if a not in "M"]
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


@dataclass
class PlantedLocus:
    locus_id: str
    contig: str
    strand: str  # "+" | "-"
    motif: str
    aa_len: int
    has_array: bool
    array_copies: int
    kind: str  # "true_art" | "no_array" | "short_upstream" | "dense_decoy" | "fragment"
    upstream_nt: int


def _dna(rng: random.Random, n: int, gc: float = 0.42) -> str:
    return "".join(rng.choices("ACGT", weights=[(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2],
                               k=n))


def _protein(rng: random.Random, n_aa: int, motif: str) -> str:
    body = [rng.choice(_AA_NO_MOTIF_CHARS) for _ in range(n_aa - 1)]
    at = n_aa // 2
    body[at:at + len(motif)] = list(motif)
    return "M" + "".join(body[: n_aa - 1])


def _orf(rng: random.Random, protein: str) -> str:
    return "".join(rng.choice(_CODON[a]) for a in protein) + rng.choice(["TAA", "TAG", "TGA"])


def _array(rng: random.Random, copies: int, unit_len: int, spacing: tuple[int, int],
           mismatch_copy: int | None = None) -> str:
    unit = _dna(rng, unit_len, gc=0.5)
    out: list[str] = []
    for i in range(copies):
        u = unit
        if mismatch_copy == i:
            j = unit_len // 2
            u = u[:j] + ("A" if u[j] != "A" else "C") + u[j + 1:]
        out.append(u)
        if i < copies - 1:
            gap = rng.randint(*spacing) - unit_len
            out.append(_dna(rng, gap))
    return "".join(out)


# (kind, motif, strand, aa_len, upstream_nt, array_copies, unit_len, spacing)
_PLAN = [
    ("true_art", "YADD", "+", 340, 2300, 7, 36, (130, 190)),
    ("true_art", "YADD", "-", 360, 2200, 5, 36, (140, 200)),
    ("true_art", "YVDD", "+", 330, 2600, 9, 36, (120, 180)),
    ("no_array", "FADD", "+", 350, 2300, 0, 0, (0, 0)),
    ("no_array", "YADD", "-", 345, 2200, 0, 0, (0, 0)),
    ("short_upstream", "YVDD", "+", 335, 700, 0, 0, (0, 0)),
    ("dense_decoy", "YADD", "+", 340, 2300, 8, 20, (36, 44)),
    ("fragment", "YADD", "+", 150, 2300, 0, 0, (0, 0)),
    ("true_art", "FADD", "-", 355, 2400, 4, 36, (150, 210)),
]


def build_world(out: Path, seed: int = 20260929) -> list[PlantedLocus]:
    """Write ``genomes.fna`` (agent-visible) and ``truth.json`` (scorer-only)."""
    rng = random.Random(seed)
    out.mkdir(parents=True, exist_ok=True)
    records: list[str] = []
    truth: list[PlantedLocus] = []
    for i, (kind, motif, strand, aa_len, up_nt, copies, unit, spacing) in enumerate(_PLAN, 1):
        contig = f"phage_{i:02d}"
        locus_id = f"L{i:02d}"
        if copies:
            arr = _array(rng, copies, unit, spacing, mismatch_copy=2 if i == 3 else None)
            lead = _dna(rng, rng.randint(150, 350))
            tail = up_nt - len(lead) - len(arr)
            upstream = lead + arr + _dna(rng, max(tail, 200))
        else:
            upstream = _dna(rng, up_nt)
        gene = _orf(rng, _protein(rng, aa_len, motif))
        block = upstream + gene
        if strand == "-":
            block = revcomp(block)
        # the short-upstream decoy sits near the contig end, so the RT truly has
        # less than 1,500 nt of sequence 5' of it (paper: "not assessed")
        head = rng.randint(20, 80) if kind == "short_upstream" else rng.randint(600, 1200)
        contig_seq = _dna(rng, head) + block + _dna(rng, rng.randint(800, 1500))
        records.append(f">{contig}\n" + "\n".join(
            contig_seq[j:j + 80] for j in range(0, len(contig_seq), 80)))
        truth.append(PlantedLocus(
            locus_id=locus_id, contig=contig, strand=strand, motif=motif, aa_len=aa_len,
            has_array=kind == "true_art", array_copies=copies if kind == "true_art" else 0,
            kind=kind, upstream_nt=len(upstream)))
    (out / "genomes.fna").write_text("\n".join(records) + "\n")
    (out / "truth.json").write_text(json.dumps([asdict(t) for t in truth], indent=2) + "\n")
    return truth


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}
