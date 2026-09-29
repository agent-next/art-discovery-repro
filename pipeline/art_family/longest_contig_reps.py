#!/usr/bin/env python3
"""Pick each MMseqs2 cluster's representative as the member on the longest contig.

paper Methods "Family definition" p.31: the 823 full-length RTs were clustered with
MMseqs2 easy-cluster into 230 clusters, "and the member on the longest contig was
taken as the representative of each cluster". MMseqs' own representative follows
similarity/length heuristics, so it is replaced here.

Usage: longest_contig_reps.py CLUSTER_TSV CONTIG_LEN_TSV OUT_IDS
  CLUSTER_TSV     mmseqs easy-cluster `<prefix>_cluster.tsv`: representative<TAB>member
  CONTIG_LEN_TSV  protein_id<TAB>length of the contig encoding it
                  (NOT-IN-PAPER: file format; the lengths come from the database
                  annotation the paper used)
  OUT_IDS         one chosen member id per cluster, ordered by MMseqs representative

NOT-IN-PAPER: ties on contig length go to the lexicographically smallest id (the
paper does not say).
"""

from __future__ import annotations

import sys
from pathlib import Path


def read_clusters(path: Path) -> dict[str, list[str]]:
    """MMseqs representative -> members (representative included), in file order."""
    clusters: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            rep, member = line.rstrip("\n").split("\t")[:2]
            clusters.setdefault(rep, []).append(member)
    return clusters


def read_contig_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            protein, length = line.rstrip("\n").split("\t")[:2]
            lengths[protein] = int(length)
    return lengths


def pick_representatives(clusters: dict[str, list[str]],
                         contig_len: dict[str, int]) -> list[str]:
    """One member per cluster: longest contig, ties to the smallest id.

    Raises KeyError listing every member without a contig length: silently ranking
    an unknown contig as 0 would change which sequence represents the cluster.
    """
    missing = sorted({m for members in clusters.values() for m in members
                      if m not in contig_len})
    if missing:
        raise KeyError(f"{len(missing)} cluster members lack a contig length: "
                       f"{', '.join(missing[:10])}")
    return [min(members, key=lambda m: (-contig_len[m], m))
            for members in clusters.values()]


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    clusters = read_clusters(Path(argv[1]))
    reps = pick_representatives(clusters, read_contig_lengths(Path(argv[2])))
    Path(argv[3]).write_text("\n".join(reps) + "\n")
    print(f"{len(reps)} clusters -> {len(reps)} representatives", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
