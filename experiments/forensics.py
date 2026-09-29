"""Replicate-campaign forensics — paper Methods "Replicate campaigns" p.38.

After rerunning the campaign, the paper asked: did any replicate rediscover the ART
loci? Procedure (parameter-exact):

- collect identifiers from the reference campaign: RT ids (paper: 130), contig ids
  (paper: 171), and names/accessions of cultured phages encoding such RTs;
- search all replicate task records (paper: 3,084 records) and session transcripts
  (paper: 5,632) for those identifiers;
- for every session that named a relevant contig (paper: 11), parse events in order
  and record (a) tool results containing a contiguous DNA string of >=200 nt and
  (b) remarks on repeats;
- an ART locus on a contig outside the identifier set would not be found by this
  search (paper's own caveat).

USAGE
    python3 experiments/forensics.py --records <campaign-root> \\
        --identifiers ids.txt --out findings.json

ids.txt: one identifier per line (RT task labels, contig names, accessions).
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

DNA_RUN = re.compile(r"[ACGTacgt]{200,}")  # paper: contiguous DNA string of >=200 nt
REPEAT_WORDS = re.compile(
    r"\b(repeats?|tandem|arrays?|direct repeats?|spacers?)\b",
    re.IGNORECASE)  # plurals: "repeats upstream of the RT" was never matched
# (S1 finding: only singular forms existed, so the headline signal undercounted)


def _unwrap_dna_lines(lines: list[str]) -> list[str]:
    """Rejoin wrapped FASTA/sequence output: pure-sequence lines are merged so a
    3,000-nt flank wrapped at 60 columns still yields one >=200-nt run (S1
    finding: line-based matching reported dna_runs=0 for exactly the sessions
    that had read the whole flank). Non-sequence lines stay as separators, so
    runs can never span across a remark or other text."""
    out: list[str] = []
    buf: list[str] = []
    for ln in lines:
        stripped = ln.strip()
        if stripped and re.fullmatch(r"[ACGTNacgtn]+", stripped):
            buf.append(stripped)
            continue
        if buf:
            out.append("".join(buf))
            buf = []
        out.append(ln)
    if buf:
        out.append("".join(buf))
    return out


def load_identifiers(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


# paper Methods p.38: 130 RT ids, 171 contig ids
PAPER_IDENTIFIER_COUNTS = {"rt_ids": 130, "contig_ids": 171}


@dataclass(frozen=True)
class IdentifierSets:
    """The paper's TWO identifier sets (Methods p.38: 130 RT ids, 171 contig
    ids) — kept distinct because only sessions naming a CONTIG are parsed
    event-by-event; a flat set cannot express that."""

    rt_ids: frozenset[str]
    contig_ids: frozenset[str]

    @property
    def paper_counts(self) -> dict[str, int]:
        return dict(PAPER_IDENTIFIER_COUNTS)  # paper p.38

    @property
    def all(self) -> frozenset[str]:
        return self.rt_ids | self.contig_ids

    def contig_named(self, named: set[str]) -> set[str]:
        """Which of the identifiers a session named are CONTIG ids."""
        return {n for n in named if n in self.contig_ids}


def load_identifier_sets(rt_path: Path, contig_path: Path) -> IdentifierSets:
    return IdentifierSets(
        rt_ids=frozenset(load_identifiers(rt_path)),
        contig_ids=frozenset(load_identifiers(contig_path)))


def scan_records(records_root: Path, ids: set[str]) -> dict[str, list[str]]:
    """task id -> identifiers mentioned anywhere in its record files."""
    hits: dict[str, list[str]] = {}
    for task_dir in sorted(Path(records_root).glob("t*")):
        if not task_dir.is_dir():
            continue
        text = "\n".join(p.read_text(errors="replace")
                         for p in task_dir.rglob("*") if p.is_file())
        found = sorted(i for i in ids if re.search(rf"\b{re.escape(i)}\b", text))
        if found:
            hits[task_dir.name] = found
    return hits


def scan_transcripts(transcripts_root: Path, ids: set[str]) -> list[dict]:
    """Per transcript naming an identifier: ordered event walk — a repeat remark
    counts as 'downstream of DNA retrieval' only when a >=200-nt DNA run appeared
    earlier in the transcript (paper p.38; grok round-2 finding 9). Identifiers
    match on word boundaries, not bare substrings."""
    out: list[dict] = []
    for tf in sorted(Path(transcripts_root).rglob("*")):
        if not tf.is_file():
            continue
        lines = _unwrap_dna_lines(
            tf.read_text(errors="replace").splitlines())
        named = sorted(i for i in ids
                       if re.search(rf"\b{re.escape(i)}\b", "\n".join(lines)))
        if not named:
            continue
        dna_seen_at: list[int] = []
        repeat_remarks = 0
        remark_after_dna = 0
        samples: list[str] = []
        for ln_no, ln in enumerate(lines):
            dna_m = DNA_RUN.search(ln)
            rem_m = REPEAT_WORDS.search(ln)
            # DNA from strictly earlier lines precedes anything on this line
            prior_dna_lines = [n for n in dna_seen_at if n < ln_no]
            if dna_m:
                dna_seen_at.append(ln_no)
            if rem_m:
                repeat_remarks += 1
                # after-DNA = DNA on an earlier line, or a DNA run starting
                # before the remark on the same line (grok round-4 finding 1:
                # the current line's own DNA must not shadow earlier lines)
                after = bool(prior_dna_lines) or (
                    dna_m is not None and dna_m.start() < rem_m.start())
                if after:
                    remark_after_dna += 1
                    if len(samples) < 3:
                        samples.append(ln.strip()[:120])
        out.append({
            "transcript": str(tf),
            "identifiers": named,
            "dna_runs_ge200nt": len(dna_seen_at),
            "repeat_remarks": repeat_remarks,
            "repeat_remarks_after_dna": remark_after_dna,
            "repeat_remark_samples": samples,
        })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--records", required=True,
                    help="campaign root containing records/tNNNN dirs")
    ap.add_argument("--transcripts", default=None,
                    help="optional root of session transcripts (rglob)")
    ap.add_argument("--identifiers", required=True, help="one identifier per line")
    ap.add_argument("--out", required=True, help="findings JSON path")
    args = ap.parse_args()

    ids = load_identifiers(Path(args.identifiers))
    records_root = Path(args.records)
    # the help says "campaign root containing records/"; older invocations that
    # pass the campaign PARENT silently produced record_hits=0 (S1 finding)
    if not records_root.exists():
        ap.error(f"--records path does not exist: {records_root}")
    if records_root.name != "records" and (records_root / "records").is_dir():
        records_root = records_root / "records"
    report = {
        "identifiers_searched": len(ids),
        "records_root": str(records_root),
        "record_hits": scan_records(records_root, ids),
        "transcript_hits": (scan_transcripts(Path(args.transcripts), ids)
                            if args.transcripts else []),
    }
    Path(args.out).write_text(json.dumps(report, indent=2) + "\n")
    print(f"identifiers={len(ids)} record_hits={len(report['record_hits'])} "
          f"transcript_hits={len(report['transcript_hits'])} -> {args.out}")


if __name__ == "__main__":
    main()
