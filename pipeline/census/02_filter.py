"""Census step 2 — filter hmmsearch domain hits to full-length RT candidates.

paper Methods "RT census by the agents" p.29: retain hits covering the complete
RT core profile (coverage >= 0.75) at a class-specific minimum length of
225-346 aa; discard fragments (n=2,016,377), weak hits (n=1,171,334; jointly bitscore
< 25, coverage < 0.35, and no YxDD catalytic motif), or excluded for other
reasons (n=622). Thresholds were fixed beforehand on the five reference RTs
named in the brief (Ec86, LtrA, BPP-1 Brt, AbiK, RT-Cas1 fusion) and on 2,339
labeled myRT sequences, retaining 99.0% of full-length labeled RTs.
Expected retained: 203,381 of 3,391,714 candidates (paper census table).

Input is an HMMER3 ``--domtblout`` whitespace table (step 01). Stdlib only.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path

# paper: expected counts around this filter (p.29)
EXPECTED_INPUT_HITS = 3_391_714
EXPECTED_DISCARDS = {"fragments": 2_016_377, "weak": 1_171_334, "other": 622}
EXPECTED_RETAINED = 203_381

# paper: "complete RT core profile" coverage floor
MIN_CORE_COVERAGE = 0.75
# paper: weak-hit thresholds (bitscore < 25 AND coverage < 0.35 AND no YxDD,
# all three jointly — Methods p.29)
WEAK_BITSCORE = 25.0
WEAK_COVERAGE = 0.35
# paper: RTs carry a YxDD active-site motif; motif-less hits are discarded
YXDD_MOTIF = re.compile(r"Y.DD")

# paper: class-specific minimum length is 225-346 aa.
# PLACEHOLDER-CLASS-SPECIFIC: the paper gives only the 225-346 aa band, not per-
# placeholders below pin every class to the band floor (225 aa) until the real
# per-class table is recovered.
CLASS_MIN_LENGTH_AA: dict[str, int] = {
    "retron": 225,
    "UG": 225,
    "group II intron": 225,
    "DGR": 225,
    "CRISPR-assoc": 225,
    "group II-like": 225,
    "Abi": 225,
    "novel": 225,
    "unknown": 225,
}
DEFAULT_MIN_LENGTH_AA = 225  # paper: lower bound of the stated 225-346 aa band

# NOT-IN-PAPER: profile-name -> class keyword map. The paper does not publish
# the 52 profile names, so class assignment from the profile name is by
# case-insensitive substring, most specific first.
CLASS_KEYWORDS: tuple[tuple[str, str], ...] = (
    ("group_ii_like", "group II-like"),
    ("group-ii-like", "group II-like"),
    ("group_ii", "group II intron"),
    ("group-ii", "group II intron"),
    ("intron", "group II intron"),
    ("crispr", "CRISPR-assoc"),
    ("retron", "retron"),
    ("dgr", "DGR"),
    ("abi", "Abi"),
    ("ug", "UG"),
    ("novel", "novel"),
)


@dataclass
class Hit:
    """One protein hit: the best domain line per target in the domtbl."""

    target: str
    tlen: int  # target (protein) length, aa
    query: str  # profile-HMM name
    qlen: int  # profile-HMM length, consensus positions
    bitscore: float  # best domain bitscore
    coverage: float  # best domain envelope coverage of the profile (0-1]
    rt_class: str = "unknown"
    # GAP: domtblout carries no sequence, so the YxDD screen needs an optional
    # FASTA of hit sequences (--seqs). None = not tested, treated as pass.
    has_yxdd: bool | None = None


@dataclass(frozen=True)
class Thresholds:
    min_core_coverage: float = MIN_CORE_COVERAGE
    weak_bitscore: float = WEAK_BITSCORE
    weak_coverage: float = WEAK_COVERAGE
    class_min_length: Mapping[str, int] = field(
        default_factory=lambda: dict(CLASS_MIN_LENGTH_AA)
    )
    default_min_length: int = DEFAULT_MIN_LENGTH_AA

    def min_length_for(self, rt_class: str) -> int:
        return self.class_min_length.get(rt_class, self.default_min_length)


def classify_profile(query_name: str) -> str:
    """Map a profile-HMM name to an RT class (NOT-IN-PAPER keyword match)."""
    name = query_name.lower()
    for needle, cls in CLASS_KEYWORDS:
        if needle in name:
            return cls
    return "unknown"


def parse_hmmsearch_domtbl(path: str | Path) -> Iterator[Hit]:
    """Yield one Hit per target: its best profile by best-domain bitscore, with
    coverage taken over the union of that profile's domain lines.

    HMMER3 domtblout columns used: target name (0), tlen (2), query name (3),
    qlen (5), domain bitscore (13), hmm from/to (15/16) — model coordinates —
    for profile coverage (envelope cols 20/21 are sequence coords; grok
    review 2026-09-24).
    """
    best_dom: dict[tuple[str, str], Hit] = {}
    spans: dict[tuple[str, str], list[tuple[int, int]]] = {}
    with Path(path).open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                continue
            qlen = int(f[5])  # HMM model length (hmmsearch query)
            hmm_from, hmm_to = int(f[15]), int(f[16])
            hit = Hit(
                target=f[0],
                tlen=int(f[2]),
                query=f[3],
                qlen=qlen,
                bitscore=float(f[13]),
                # profile coverage = aligned model columns / model length
                # (grok review 2026-09-24: envelope cols 20/21 are sequence
                # coordinates, not profile coordinates)
                coverage=(hmm_to - hmm_from + 1) / qlen,
                rt_class=classify_profile(f[3]),
            )
            key = (hit.target, hit.query)
            spans.setdefault(key, []).append((hmm_from, hmm_to))
            prev = best_dom.get(key)
            if prev is None or hit.bitscore > prev.bitscore:
                best_dom[key] = hit
    best: dict[str, Hit] = {}
    for key, hit in best_dom.items():
        # An RT core split by an insertion is reported as several domains of one
        # profile; the core is covered when their model columns jointly are.
        hit.coverage = _union_length(spans[key]) / hit.qlen
        prev = best.get(hit.target)
        if prev is None or hit.bitscore > prev.bitscore:
            best[hit.target] = hit
    yield from best.values()


def _union_length(intervals: list[tuple[int, int]]) -> int:
    total, end = 0, 0
    for lo, hi in sorted(intervals):
        lo = max(lo, end + 1)
        if hi >= lo:
            total += hi - lo + 1
            end = hi
    return total


def is_weak_hit(hit: Hit, thresholds: Thresholds) -> bool:
    """paper (Methods p.29): weak hits = bit score < 25 AND coverage < 0.35 AND
    no YxDD catalytic motif — all three simultaneously (grok review 2026-09-24:
    the OR reading discarded candidates the paper keeps)."""
    if (hit.bitscore < thresholds.weak_bitscore
            and hit.coverage < thresholds.weak_coverage
            and hit.has_yxdd is False):
        return True
    return False


def passes_filters(hit: Hit, thresholds: Thresholds | None = None) -> bool:
    """True when the hit survives every step-02 census filter."""
    th = thresholds or Thresholds()
    if is_weak_hit(hit, th):  # paper: discard weak hits
        return False
    if hit.coverage < th.min_core_coverage:  # paper: discard fragments (<0.75 core cov)
        return False
    return hit.tlen >= th.min_length_for(hit.rt_class)  # paper: class min length


def load_yxdd_fasta(path: str | Path) -> dict[str, bool]:
    """Map sequence id -> has-YxDD-motif from a FASTA of candidate proteins."""
    flags: dict[str, bool] = {}
    name: str | None = None
    chunks: list[str] = []
    with Path(path).open() as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    flags[name] = bool(YXDD_MOTIF.search("".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        flags[name] = bool(YXDD_MOTIF.search("".join(chunks)))
    return flags


def plan() -> str:
    lines = [
        "census step 02 -- filter hmmsearch hits (dry-run plan)",
        f"  input: domtblout from step 01 (expected ~{EXPECTED_INPUT_HITS:,} candidates)",
        f"  retain: core profile coverage >= {MIN_CORE_COVERAGE}"
        " AND length >= class minimum (paper band 225-346 aa)",
        "  class min lengths: "
        + ", ".join(f"{k}={v}" for k, v in sorted(CLASS_MIN_LENGTH_AA.items())),
        "  GAP: per-class minima not stated -- placeholders at band floor 225",
        f"  discard weak: bitscore < {WEAK_BITSCORE} AND coverage < {WEAK_COVERAGE}"
        " AND no YxDD motif, jointly (needs --seqs FASTA to test)",
        f"  expected retained: {EXPECTED_RETAINED:,}",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("domtbl", nargs="?", help="HMMER3 --domtblout file from step 01")
    ap.add_argument("--seqs", help="optional FASTA of hit sequences for the YxDD check")
    ap.add_argument("--out", help="retained target ids (default: stdout)")
    ap.add_argument("--dry-run", action="store_true", help="print the step plan and exit")
    args = ap.parse_args(argv)

    if args.dry_run:
        print(plan())
        return 0
    if not args.domtbl:
        ap.error("domtbl is required unless --dry-run")

    yxdd = load_yxdd_fasta(args.seqs) if args.seqs else {}
    th = Thresholds()
    kept: list[str] = []
    n_seen = 0
    for hit in parse_hmmsearch_domtbl(args.domtbl):
        n_seen += 1
        if hit.target in yxdd:
            hit.has_yxdd = yxdd[hit.target]
        if passes_filters(hit, th):
            kept.append(hit.target)

    out = Path(args.out).open("w") if args.out else sys.stdout
    try:
        for target in kept:
            out.write(target + "\n")
    finally:
        if args.out:
            out.close()
    print(
        f"retained {len(kept)} of {n_seen} hits (paper: {EXPECTED_RETAINED:,})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

