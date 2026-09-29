"""Build the RNA-seq feature file from a GenBank record.

paper Methods p.36-37: features = 258 annotated SA1 CDS + the array RNA (259 total).
NOT-IN-PAPER: the paper does not publish its feature file, so it is rebuilt here from
the MW218148.1 GenBank CDS annotation; the array RNA is added as one extra feature
per strand the caller supplies (for SA1, the locus called by
``scripts/scan_genome.py``). The paper does not say which strand carries the array RNA,
so passing both strands reports both; the 259-feature set uses the stronger one. Features are typed ``CDS`` so ``featureCounts -t CDS``
counts them. Stdlib only.

Usage: genbank_to_gtf.py IN.gb OUT.gtf [--array START END STRAND]...
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

_LOC = re.compile(r"^(complement\()?[<>]?(\d+)\.\.[<>]?(\d+)\)?$")


def parse_cds(text: str) -> list[tuple[str, int, int, str]]:
    """(id, start, end, strand) per CDS, 1-based inclusive. Joined locations are refused."""
    feats: list[tuple[str, int, int, str]] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("     CDS "):
            loc = line[21:].strip()
            i += 1
            while (i < len(lines) and lines[i].startswith(" " * 21)
                   and not lines[i].strip().startswith("/")):
                loc += lines[i].strip()
                i += 1
            m = _LOC.match(loc)
            if not m:
                raise ValueError(f"unsupported CDS location: {loc}")
            ident = None
            while i < len(lines) and lines[i].startswith(" " * 21):
                q = re.match(r'\s+/(protein_id|locus_tag)="([^"]+)"', lines[i])
                if q and (ident is None or q.group(1) == "protein_id"):
                    ident = q.group(2)
                i += 1
            feats.append((ident or f"cds{len(feats) + 1}", int(m.group(2)), int(m.group(3)),
                          "-" if m.group(1) else "+"))
            continue
        i += 1
    return feats


def seqid(text: str) -> str:
    m = re.search(r"^VERSION\s+(\S+)", text, re.M)
    if not m:
        raise ValueError("no VERSION line")
    return m.group(1)


def to_gtf(text: str, arrays: list[tuple[int, int, str]] | None = None) -> str:
    sid = seqid(text)
    rows = parse_cds(text)
    for start, end, strand in arrays or []:
        rows.append((f"ART_array_RNA_{'fwd' if strand == '+' else 'rev'}", start, end, strand))
    ids = [r[0] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate feature ids")
    return "".join(
        f'{sid}\tart-harness\tCDS\t{s}\t{e}\t.\t{st}\t0\tgene_id "{fid}"; transcript_id "{fid}";\n'
        for fid, s, e, st in rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("genbank")
    ap.add_argument("out")
    ap.add_argument("--array", nargs=3, action="append", metavar=("START", "END", "STRAND"))
    args = ap.parse_args(argv)
    arrays = []
    for start, end, strand in args.array or []:
        if strand not in ("+", "-"):
            ap.error("STRAND must be + or -")
        arrays.append((int(start), int(end), strand))
    gtf = to_gtf(Path(args.genbank).read_text(), arrays)
    Path(args.out).write_text(gtf)
    print(f"{gtf.count(chr(10))} features -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
