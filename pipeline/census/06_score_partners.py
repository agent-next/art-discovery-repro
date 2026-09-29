"""Census step 6 — partner-gene family scoring and promotion.

paper Methods "Scoring of partner families" p.30: neighboring proteins were
grouped into families by best Pfam-A match, else MMseqs2 clustering; the 3,564
families recurring across loci were evaluated against three promotion filters
(thresholds set by the census worker before the census ran):
  (1) conserved within >=1 RT clade (>=10% of the clade's loci OR >=3 loci);
  (2) independent occurrences: >=3 clusters at 90% identity, >=3 biosamples,
      beside RTs of >=2 classes;
  (3) genes lie nearer to the RT than randomly drawn genes of the same loci
      (permutation P <= 0.05) AND lie on the strand of the RT OR within 100 bp
      of the adjacent gene.
paper: known RT partners and general mobile-element residents were designated
controls beforehand; 46 families met all three filters (13 of them controls
incl. Avd, Cas1, Cas2), 16 were promoted after annotation review, and a 17th
via a follow-up task.

GAP: the paper does not publish the permutation null model or the proximity
statistic; here the statistic is the count of occurrences passing
proximity_qualifies() against a draw from the neighborhood-gene pool.

Pure functions over simple record types; the permutation test is seeded.
Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import random
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

# paper: promotion filter 1 -- conserved within >=1 RT clade
CONSERVED_MIN_FRACTION = 0.10  # >=10% of the clade's loci
CONSERVED_MIN_LOCI = 3  # OR >=3 loci

# paper: promotion filter 2 -- independent occurrences
INDEPENDENT_MIN_CLUSTERS_90 = 3  # >=3 clusters at 90% identity
INDEPENDENT_MIN_BIOSAMPLES = 3  # >=3 biosamples
INDEPENDENT_MIN_RT_CLASSES = 2  # beside RTs of >=2 classes

# paper: promotion filter 3 -- permutation test P <= 0.05 for RT proximity
# AND (same strand OR within 100 bp of the adjacent gene)
PERMUTATION_P_MAX = 0.05
ADJACENT_MAX_BP = 100
# GAP: the paper does not state the permutation count; 10,000 placeholder.
N_PERMUTATIONS = 10_000
# NOT-IN-PAPER: the paper states no seed; fixed for determinism.
DEFAULT_SEED = 1

# paper: control families were designated beforehand.
# GAP: the paper does not list them; populate before running for real.
CONTROL_FAMILIES: frozenset[str] = frozenset()


@dataclass(frozen=True)
class NeighborGene:
    """One occurrence of a candidate partner family in an RT neighborhood."""

    gene_id: str
    locus_id: str
    family: str  # from assign_family()
    rt_clade: str  # RT clade of the neighborhood's anchor
    rt_class: str  # class of the neighboring RT
    cluster_id_90: str  # protein cluster at 90% identity
    biosample: str
    rt_adjacent: bool  # directly next to the RT gene
    same_strand_as_rt: bool
    distance_to_adjacent_bp: int | None  # to the adjacent gene on the RT side


def assign_family(best_pfam_family: str | None, mmseqs_cluster: str) -> str:
    """paper: family assignment by best Pfam-A match, else MMseqs2 cluster."""
    return best_pfam_family if best_pfam_family else f"cluster:{mmseqs_cluster}"


def conserved_in_clade(loci_with_family: int, clade_loci_total: int) -> bool:
    """paper: >=10% of the clade's loci OR >=3 loci."""
    return loci_with_family >= CONSERVED_MIN_LOCI or (
        clade_loci_total > 0
        and loci_with_family / clade_loci_total >= CONSERVED_MIN_FRACTION
    )


def passes_conservation(
    occurrences: Sequence[NeighborGene], clade_loci_totals: Mapping[str, int]
) -> bool:
    """Filter 1: conserved within >=1 RT clade."""
    by_clade: dict[str, set[str]] = defaultdict(set)
    for occ in occurrences:
        by_clade[occ.rt_clade].add(occ.locus_id)
    return any(
        conserved_in_clade(len(loci), clade_loci_totals.get(clade, 0))
        for clade, loci in by_clade.items()
    )


def passes_independence(occurrences: Sequence[NeighborGene]) -> bool:
    """Filter 2: >=3 90%-id clusters, >=3 biosamples, >=2 RT classes."""
    return (
        len({o.cluster_id_90 for o in occurrences}) >= INDEPENDENT_MIN_CLUSTERS_90
        and len({o.biosample for o in occurrences}) >= INDEPENDENT_MIN_BIOSAMPLES
        and len({o.rt_class for o in occurrences}) >= INDEPENDENT_MIN_RT_CLASSES
    )


def proximity_qualifies(occ: NeighborGene) -> bool:
    """paper: RT proximity AND (same strand OR within 100 bp of adjacent gene)."""
    if not occ.rt_adjacent:
        return False
    return occ.same_strand_as_rt or (
        occ.distance_to_adjacent_bp is not None
        and occ.distance_to_adjacent_bp <= ADJACENT_MAX_BP
    )


def permutation_pvalue(
    occurrences: Sequence[NeighborGene],
    background: Sequence[NeighborGene],
    n_permutations: int = N_PERMUTATIONS,
    rng: random.Random | None = None,
) -> float:
    """Seeded permutation test on the qualifying-proximity count.

    Paper p.30 (verbatim): "genes of the family had to lie nearer to the RT
    than randomly drawn genes of the same loci (P <= 0.05)". The NULL POOL is
    stated (genes of the same loci — the caller restricts `background`); the
    DISTANCE MEASURE is not, and graded distance-to-RT is not in the neighbor
    records, so the STATISTIC is INTERPRETED as the count of occurrences
    passing proximity_qualifies() (RT-adjacent AND same-strand-or-<=100bp).
    Draw mechanics (k draws without replacement) also unstated. GAP: exact
    nearness statistic — needs a distance-to-RT field from the producer.
    """
    rng = rng or random.Random(DEFAULT_SEED)
    observed = sum(1 for o in occurrences if proximity_qualifies(o))
    k = len(occurrences)
    pool = list(background)
    exceed = 0
    for _ in range(n_permutations):
        draw = rng.sample(pool, min(k, len(pool)))
        if sum(1 for o in draw if proximity_qualifies(o)) >= observed:
            exceed += 1
    return (exceed + 1) / (n_permutations + 1)


def promotes(
    occurrences: Sequence[NeighborGene],
    clade_loci_totals: Mapping[str, int],
    background: Sequence[NeighborGene],
    n_permutations: int = N_PERMUTATIONS,
    rng: random.Random | None = None,
) -> tuple[bool, float]:
    """A family is promoted when all three paper filters pass.

    Returns (promoted, permutation_pvalue).
    """
    if not passes_conservation(occurrences, clade_loci_totals):
        return False, 1.0
    if not passes_independence(occurrences):
        return False, 1.0
    # Paper p.30: the null draws come from "genes of the same loci" — restrict
    # the background to the family's loci before the permutation test (the
    # full neighborhood pool made distant loci do the null's job).
    family_loci = {o.locus_id for o in occurrences}
    same_loci_pool = [g for g in background if g.locus_id in family_loci]
    pool = same_loci_pool if same_loci_pool else list(background)
    p = permutation_pvalue(occurrences, pool, n_permutations, rng)
    return p <= PERMUTATION_P_MAX, p


def plan() -> str:
    lines = [
        "census step 06 -- partner-gene family scoring (dry-run plan)",
        "  family: best Pfam-A match, else MMseqs2 cluster",
        "  promotion filters (paper):",
        f"    1. conserved in >=1 RT clade: >={CONSERVED_MIN_FRACTION:.0%} of loci"
        f" OR >= {CONSERVED_MIN_LOCI} loci",
        f"    2. independent occurrences: >= {INDEPENDENT_MIN_CLUSTERS_90} clusters"
        f" at 90% id, >= {INDEPENDENT_MIN_BIOSAMPLES} biosamples, RTs of"
        f" >= {INDEPENDENT_MIN_RT_CLASSES} classes",
        f"    3. permutation test P <= {PERMUTATION_P_MAX} for RT proximity AND"
        f" (same strand OR within {ADJACENT_MAX_BP} bp of adjacent gene);"
        " null pool: genes of the SAME LOCI (paper p.30)",
        f" n={N_PERMUTATIONS} permutations (GAP placeholder)",
        "  controls: designated beforehand (GAP: identities not listed in paper)",
    ]
    return "\n".join(lines)


def load_neighbors(path: str | Path) -> list[NeighborGene]:
    """Read TSV: gene_id, locus_id, family, rt_clade, rt_class, cluster_id_90,
    biosample, rt_adjacent(0/1), same_strand_as_rt(0/1), distance_to_adjacent_bp(-).

    NOT-IN-PAPER: TSV plumbing; the paper fixes the filters, not file formats.
    """
    out: list[NeighborGene] = []
    with Path(path).open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#") or len(row) < 10:
                continue
            out.append(
                NeighborGene(
                    gene_id=row[0],
                    locus_id=row[1],
                    family=row[2],
                    rt_clade=row[3],
                    rt_class=row[4],
                    cluster_id_90=row[5],
                    biosample=row[6],
                    rt_adjacent=row[7] == "1",
                    same_strand_as_rt=row[8] == "1",
                    distance_to_adjacent_bp=None if row[9] in ("-", "") else int(row[9]),
                )
            )
    return out


def load_clade_totals(path: str | Path) -> dict[str, int]:
    """Read TSV: rt_clade, loci_total. NOT-IN-PAPER plumbing."""
    totals: dict[str, int] = {}
    with Path(path).open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#") or len(row) < 2:
                continue
            totals[row[0]] = int(row[1])
    return totals


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("neighbors", nargs="?", help="NeighborGene TSV (see load_neighbors)")
    ap.add_argument("--clade-totals", help="TSV: rt_clade, loci_total")
    ap.add_argument("--permutations", type=int, default=N_PERMUTATIONS,
                    help=f"permutation count (default {N_PERMUTATIONS}; GAP placeholder)")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED,
                    help=f"permutation seed (default {DEFAULT_SEED}; NOT-IN-PAPER)")
    ap.add_argument("--dry-run", action="store_true", help="print the step plan and exit")
    args = ap.parse_args(argv)

    if args.dry_run:
        print(plan())
        return 0
    if not args.neighbors:
        ap.error("neighbors TSV is required unless --dry-run")
    if not args.clade_totals:
        ap.error("--clade-totals is required unless --dry-run")

    genes = load_neighbors(args.neighbors)
    totals = load_clade_totals(args.clade_totals)
    by_family: dict[str, list[NeighborGene]] = defaultdict(list)
    for gene in genes:
        by_family[gene.family].append(gene)

    rng = random.Random(args.seed)
    for family in sorted(by_family):
        promoted, p = promotes(
            by_family[family], totals, genes, args.permutations, rng
        )
        control = "\tcontrol" if family in CONTROL_FAMILIES else ""
        verdict = "promoted" if promoted else "rejected"
        print(f"{family}\t{verdict}\tp={p:.4f}{control}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
