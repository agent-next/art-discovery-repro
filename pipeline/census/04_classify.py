"""Census step 4 — three-tier RT classification.

paper Methods "RT census by the agents" p.29:
tier 1: 38 lineage-specific profile HMMs built from myRT (mafft --auto,
hmmbuild); a target is assigned when its best profile clears a class-specific
bitscore floor with a margin of at least 10 bits over every other class
(precision 0.997 in five-fold cross-validation).
tier 2: MMseqs2 search against labeled RTs; accept at bitscore >= 150,
query coverage >= 0.5, identity >= 30%.
tier 3: FastTree -lg -gamma of 2,241 RT domains aligned to RVT_1; assign when
the five nearest labeled leaves agree.
paper: the tiers placed 35,168 clusters in the seven known classes;
EXPECTED_CLASS_COUNTS holds the census table.

Tier-1 margin and tier-2 threshold logic are pure functions; the FastTree tier
shells out (--dry-run prints the command). Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

# paper: 38 lineage-specific profile HMMs in tier 1
TIER1_PROFILE_COUNT = 38
# paper: tier-1 assignment needs a >=10-bit margin over the next-best profile
TIER1_MIN_MARGIN_BITS = 10.0

# GAP: the paper states class-specific bitscore floors for tier 1 but does not
# PLACEHOLDER-CLASS-SPECIFIC (paper gives no per-class floors); uniform 100 bits
# until calibrated against the paper's myRT cross-validation - do NOT treat as
# recovered.
TIER1_BITSCORE_FLOORS: dict[str, float] = {
    "retron": 100.0,
    "UG": 100.0,
    "group II intron": 100.0,
    "DGR": 100.0,
    "CRISPR-assoc": 100.0,
    "group II-like": 100.0,
    "Abi": 100.0,
    "novel": 100.0,
}
TIER1_DEFAULT_FLOOR = 100.0  # GAP placeholder, see above

# paper: tier-2 acceptance thresholds
TIER2_MIN_BITSCORE = 150.0
TIER2_MIN_COVERAGE = 0.5
TIER2_MIN_IDENTITY = 0.30  # 30%

# paper: tier-3 assignment when the 5 nearest labeled leaves agree
TIER3_NEAREST_LEAVES = 5

# paper: expected class counts after classification (census table, p.29);
# the seven known classes total 35,168 placed clusters
EXPECTED_CLASS_COUNTS: dict[str, int] = {
    "retron": 11_517,
    "UG": 10_759,
    "group II intron": 5_461,
    "DGR": 3_857,
    "CRISPR-assoc": 2_067,
    "group II-like": 827,
    "Abi": 680,
    "novel": 25_737,
    "unplaced": 137_385,
}


@dataclass(frozen=True)
class Tier1Hit:
    """One lineage-specific profile hit (tier 1)."""

    target: str
    profile: str
    rt_class: str
    bitscore: float


@dataclass(frozen=True)
class Tier2Hit:
    """Best MMseqs2 hit of one target against the labeled RT database (tier 2)."""

    target: str
    label: str
    bitscore: float
    coverage: float
    identity: float  # fraction, 0-1


def tier1_assign(
    hits: Sequence[Tier1Hit],
    floors: Mapping[str, float] | None = None,
    margin: float = TIER1_MIN_MARGIN_BITS,
) -> str | None:
    """Tier-1 class for one target's profile hits, or None.

    paper: assign when the best profile clears the class-specific bitscore
    floor with a >=10-bit margin over the next-best profile of a different
    class.
    """
    fl = floors if floors is not None else TIER1_BITSCORE_FLOORS
    if not hits:
        return None
    ordered = sorted(hits, key=lambda h: h.bitscore, reverse=True)
    best = ordered[0]
    if best.bitscore < fl.get(best.rt_class, TIER1_DEFAULT_FLOOR):
        return None
    next_diff = next((h for h in ordered[1:] if h.rt_class != best.rt_class), None)
    if next_diff is not None and best.bitscore - next_diff.bitscore < margin:
        return None
    return best.rt_class


def tier2_accept(hit: Tier2Hit) -> bool:
    """paper: accept at bitscore >= 150, coverage >= 0.5, identity >= 30%."""
    return (
        hit.bitscore >= TIER2_MIN_BITSCORE
        and hit.coverage >= TIER2_MIN_COVERAGE
        and hit.identity >= TIER2_MIN_IDENTITY
    )


def tier3_assign(neighbor_classes: Sequence[str]) -> str | None:
    """paper: assign when the 5 nearest labeled leaves agree on a class."""
    nearest = list(neighbor_classes[:TIER3_NEAREST_LEAVES])
    if len(nearest) < TIER3_NEAREST_LEAVES:
        return None
    return nearest[0] if all(c == nearest[0] for c in nearest) else None


def fasttree_cmd(aln_path: str | Path, tree_out: str | Path) -> str:
    """paper: FastTree -lg -gamma placement tree (shells out; stdout -> tree)."""
    return f"FastTree -lg -gamma {aln_path} > {tree_out}"


def classify_target(
    target: str,
    tier1_hits: Sequence[Tier1Hit],
    tier2_hit: Tier2Hit | None = None,
    tier3_neighbors: Sequence[str] | None = None,
) -> tuple[str, str]:
    """Return (tier, rt_class) for one target across the three tiers."""
    assigned = tier1_assign(tier1_hits)
    if assigned is not None:
        return "tier1", assigned
    if tier2_hit is not None and tier2_accept(tier2_hit):
        return "tier2", tier2_hit.label
    if tier3_neighbors is not None:
        assigned = tier3_assign(tier3_neighbors)
        if assigned is not None:
            return "tier3", assigned
        # paper p.29: clusters that "formed a loosely organized group on the
        # tree apart from any labeled member" are the eighth class "novel";
        # "scattered clusters that could not be placed" are the ninth
        # "unplaced". INTERPRETED criterion (paper gives no rule): a majority
        # of unlabeled neighbors among the k nearest leaves = grouped apart
        # from labeled members -> novel; labeled-but-disagreeing neighborhood
        # = scattered -> unplaced. NOT-IN-PAPER: the majority threshold.
        nearest = tier3_neighbors[:TIER3_NEAREST_LEAVES]
        unlabeled = sum(1 for c in nearest if not c or c in ("unplaced", "-"))
        if unlabeled * 2 > len(nearest):
            return "tier3", "novel"
    return "unplaced", "unplaced"


def load_tier3_neighbors(path: str | Path) -> dict[str, list[str]]:
    """Read tier-3 TSV: target, neighbor_class (one row per placed neighbor).

    NOT-IN-PAPER: TSV plumbing; the neighbor classes come from the tier-1/2
    labels of sequences placed next to the unplaced target in the FastTree.
    """
    by_target: dict[str, list[str]] = defaultdict(list)
    with Path(path).open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            by_target[row[0]].extend(row[1:])
    return dict(by_target)


def load_tier1_hits(path: str | Path) -> dict[str, list[Tier1Hit]]:
    """Read tier-1 TSV: target, profile, rt_class, bitscore ('#' comments ok).

    NOT-IN-PAPER: TSV plumbing; the paper fixes thresholds, not file formats.
    """
    by_target: dict[str, list[Tier1Hit]] = defaultdict(list)
    with Path(path).open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#") or len(row) < 4:
                continue
            hit = Tier1Hit(target=row[0], profile=row[1], rt_class=row[2],
                           bitscore=float(row[3]))
            by_target[hit.target].append(hit)
    return dict(by_target)


def load_tier2_hits(path: str | Path) -> dict[str, Tier2Hit]:
    """Read tier-2 TSV: target, label, bitscore, coverage, identity.

    NOT-IN-PAPER: precomputed best-hit-per-target table (mmseqs convertalis
    plumbing); the paper fixes thresholds, not file formats.
    """
    best: dict[str, Tier2Hit] = {}
    with Path(path).open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#") or len(row) < 5:
                continue
            hit = Tier2Hit(
                target=row[0],
                label=row[1],
                bitscore=float(row[2]),
                coverage=float(row[3]),
                identity=float(row[4]),
            )
            prev = best.get(hit.target)
            if prev is None or hit.bitscore > prev.bitscore:
                best[hit.target] = hit
    return best


def plan(aln: str = "results/census/unplaced.aln.faa",
         tree: str = "results/census/unplaced.nwk") -> str:
    lines = [
        "census step 04 -- three-tier RT classification (dry-run plan)",
        f"  tier 1: {TIER1_PROFILE_COUNT} lineage-specific profile HMMs; assign when the",
        "          best profile clears its class bitscore floor with a "
        f">={TIER1_MIN_MARGIN_BITS:g}-bit margin over the next-best class",
        "          GAP: per-class floors not stated -- uniform placeholder "
        f"{TIER1_DEFAULT_FLOOR:g} bits",
        f"  tier 2: MMseqs2 vs labeled RTs; accept at bitscore >= {TIER2_MIN_BITSCORE:g},",
        f"          coverage >= {TIER2_MIN_COVERAGE:g}, identity >= {TIER2_MIN_IDENTITY:.0%}",
        f"  tier 3: FastTree -lg -gamma; assign when the {TIER3_NEAREST_LEAVES} nearest",
        "          labeled leaves agree:",
        f"          {fasttree_cmd(aln, tree)}",
        "  expected class counts:",
    ]
    lines += [f"    {k}: {v:,}" for k, v in EXPECTED_CLASS_COUNTS.items()]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tier1", help="TSV: target, profile, rt_class, bitscore")
    ap.add_argument("--tier2", help="TSV: target, label, bitscore, coverage, identity")
    ap.add_argument("--tier3", help="TSV: target, neighbor_class[, ...] from the "
                                   "unplaced FastTree placement; without it "
                                   "tier-1/tier-2 leftovers stay unplaced")
    ap.add_argument("--out", help="assignments TSV: target, tier, rt_class (default stdout)")
    ap.add_argument("--aln", default="results/census/unplaced.aln.faa",
                    help="alignment for the tier-3 FastTree placement")
    ap.add_argument("--tree", default="results/census/unplaced.nwk",
                    help="tier-3 FastTree output tree")
    ap.add_argument("--dry-run", action="store_true", help="print the step plan and exit")
    args = ap.parse_args(argv)

    if args.dry_run:
        print(plan(args.aln, args.tree))
        return 0
    if not args.tier1:
        ap.error("--tier1 is required unless --dry-run")

    t1 = load_tier1_hits(args.tier1)
    t2 = load_tier2_hits(args.tier2) if args.tier2 else {}
    # S1 finding: tier3_assign existed but was unreachable -- classify_target
    # never received tier-3 neighbors, so every tier-3 candidate landed in
    # "unplaced" and the paper's 137,385 unplaced count was not reproducible.
    t3 = load_tier3_neighbors(args.tier3) if args.tier3 else None
    targets = sorted(set(t1) | set(t2)
                     | (set(t3) if t3 is not None else set()))

    out = Path(args.out).open("w") if args.out else sys.stdout
    counts: dict[str, int] = defaultdict(int)
    try:
        for target in targets:
            tier, rt_class = classify_target(
                target, t1.get(target, []), t2.get(target),
                t3.get(target, []) if t3 is not None else None)
            counts[rt_class] += 1
            out.write(f"{target}\t{tier}\t{rt_class}\n")
    finally:
        if args.out:
            out.close()

    for cls, n in sorted(counts.items()):
        exp = EXPECTED_CLASS_COUNTS.get(cls)
        suffix = f" (paper: {exp:,})" if exp is not None else ""
        print(f"{cls}\t{n}{suffix}", file=sys.stderr)
    print(f"tier-3 placement command: {fasttree_cmd(args.aln, args.tree)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
