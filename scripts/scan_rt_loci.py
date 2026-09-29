"""Per-locus array scan of annotated RT genes in real GenBank records.

paper Methods p.32: the k-mer scan and the delimitation run on the window upstream of
each RT (3,000 nt for the scan, 6,000 nt for delimitation). This driver takes RT protein
accessions (for example the ids retained by ``pipeline/census/02_filter.py``), finds each
gene in a GenBank flat file, cuts the upstream window in the RT's orientation and calls
``artharness.arrays.kmer_scan`` and ``delimit_array``. The annotated CDS in the window
are passed as ``gene_spans`` so both p.32 delimitation rules (coding-repeat exclusion and
the >= 300-nt RT-adjacency rule) are live here (REPRODUCTION.md GAP-9).

USAGE: python3 scripts/scan_rt_loci.py RT_IDS.txt genomes.gbff[.gz] [-o OUT.tsv] [-j JOBS]
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import random
import re
import sys
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from artharness.arrays import MAX_UPSTREAM_DELIMIT, delimit_array, kmer_scan, revcomp

SEED = 20260923  # same fixed seed as scripts/scan_genome.py
_PID = re.compile(r'/protein_id="([^"]+)"')
_INT = re.compile(r"\d+")


@dataclass
class Cds:
    protein_id: str | None
    lo: int  # 1-based inclusive
    hi: int
    strand: str


@dataclass
class Record:
    accession: str
    organism: str
    cds: list[Cds]
    seq: str


def records(path: Path) -> Iterator[str]:
    """GenBank flat-file records as text, streamed."""
    opener = gzip.open if path.suffix == ".gz" else open
    buf: list[str] = []
    with opener(path, "rt") as fh:
        for line in fh:
            buf.append(line)
            if line.startswith("//"):
                yield "".join(buf)
                buf = []
    if buf and "".join(buf).strip():
        yield "".join(buf)


def parse_record(text: str) -> Record:
    acc = re.search(r"^VERSION\s+(\S+)", text, re.M)
    org = re.search(r"^  ORGANISM\s+(.+)$", text, re.M)
    features, _, origin = text.partition("\nORIGIN")
    cds: list[Cds] = []
    lines = features.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("     CDS "):
            loc = lines[i][21:].strip()
            i += 1
            while (i < len(lines) and lines[i].startswith(" " * 21)
                   and not lines[i].strip().startswith("/")):
                loc += lines[i].strip()
                i += 1
            pid = None
            while i < len(lines) and lines[i].startswith(" " * 21):
                m = _PID.search(lines[i])
                if m:
                    pid = m.group(1)
                i += 1
            nums = [int(x) for x in _INT.findall(loc)]
            if nums:
                cds.append(Cds(pid, min(nums), max(nums),
                               "-" if loc.startswith("complement") else "+"))
            continue
        i += 1
    seq = "".join(c for c in origin.split("\n//")[0] if c.isalpha()).upper()
    return Record(acc.group(1) if acc else "?", org.group(1).strip() if org else "?", cds, seq)


def upstream_window(rec: Record, rt: Cds, length: int = MAX_UPSTREAM_DELIMIT
                    ) -> tuple[str, list[tuple[int, int]]]:
    """Upstream window oriented so the RT follows it, plus the other CDS as
    (start, end) spans in window coordinates (0-based, end exclusive)."""
    if rt.strand == "+":
        a, b = max(0, rt.lo - 1 - length), rt.lo - 1
        window = rec.seq[a:b]
    else:
        a, b = rt.hi, min(len(rec.seq), rt.hi + length)
        window = revcomp(rec.seq[a:b])
    spans = []
    for c in rec.cds:
        if c is rt:
            continue
        s0, e1 = c.lo - 1, c.hi
        if e1 <= a or s0 >= b:
            continue
        s0, e1 = max(s0, a), min(e1, b)
        spans.append((s0 - a, e1 - a) if rt.strand == "+" else (b - e1, b - s0))
    return window, spans


def scan_locus(rec: Record, rt: Cds) -> dict[str, object]:
    window, spans = upstream_window(rec, rt)
    locus = rt.protein_id or "?"
    rng = random.Random(f"{SEED}:{locus}")
    call = kmer_scan(locus, window, rng)
    row: dict[str, object] = {
        "rt": locus, "record": rec.accession, "organism": rec.organism, "strand": rt.strand,
        "rt_lo": rt.lo, "rt_hi": rt.hi, "upstream_nt": len(window), "status": call.status,
        "R": call.R, "seed": call.seed or "", "delim_copies": 0, "delim_repeat": "",
        "coding_repeat": "", "rt_adjacent": ""}
    arr = delimit_array(locus, window, random.Random(f"{SEED}:d:{locus}"), gene_spans=spans)
    if arr is not None:
        row.update(delim_copies=len(arr.copy_starts), delim_repeat=arr.repeat,
                   coding_repeat=arr.coding_repeat, rt_adjacent=arr.rt_adjacent)
    return row


def _scan_record(text: str, wanted: frozenset[str]) -> list[dict[str, object]]:
    rec = parse_record(text)
    return [scan_locus(rec, c) for c in rec.cds if c.protein_id in wanted]


COLUMNS = ["rt", "record", "organism", "strand", "rt_lo", "rt_hi", "upstream_nt", "status",
           "R", "seed", "delim_copies", "delim_repeat", "coding_repeat", "rt_adjacent"]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("rt_ids")
    ap.add_argument("gbff", type=Path)
    ap.add_argument("-o", "--out", type=Path)
    ap.add_argument("-j", "--jobs", type=int, default=1)
    args = ap.parse_args(argv)
    wanted = frozenset(x.split("|")[0].strip() for x in Path(args.rt_ids).read_text().split()
                       if x.strip())
    texts = (t for t in records(args.gbff) if wanted.intersection(_PID.findall(t)))
    rows: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for got in pool.map(_scan_record, texts, itertools.repeat(wanted)):
            rows.extend(got)
    found = {str(r["rt"]) for r in rows}
    out = ["\t".join(COLUMNS)] + ["\t".join(str(r[c]) for c in COLUMNS) for r in rows]
    (args.out.write_text("\n".join(out) + "\n") if args.out else print("\n".join(out)))
    print(f"{len(rows)} loci scanned, {len(wanted - found)} ids not found in the flat file, "
          f"{sum(r['status'] == 'array' for r in rows)} with an array call", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
