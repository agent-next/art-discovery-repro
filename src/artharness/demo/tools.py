"""The sandbox tools an agent may call.

In the paper the agents had a 60-core sandbox, bioinformatics tools, a metagenomic
database, literature search and a GPU queue (Methods p.28). Here they have five small
tools over the toy world. Two of them run the paper's own array algorithms unchanged
(``artharness.arrays.kmer_scan`` = step 1, ``delimit_array`` = step 2), so the demo's
science is real even though the input is tiny.

Every tool returns a :class:`ToolResult` whose ``headline`` is one machine-checkable
line (``key=value ...``). Workers must quote it in their summary and supervisors check
that they did; this is the demo's version of "supervisors checked quantitative claims
against the artifacts" (paper Results p.2).
"""

from __future__ import annotations

import functools
import random
import re
import statistics
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from ..arrays import MAX_UPSTREAM_SCAN, delimit_array, kmer_scan
from .world import read_fasta, revcomp

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(
    (a, b, c) for a in _BASES for b in _BASES for c in _BASES)}

RT_MOTIF = re.compile(r"[YF][ADV]DD")  # NOT-IN-PAPER: stand-in for the 52 RT HMMs
FRAGMENT_MAX_AA = 250  # NOT-IN-PAPER: toy full-length cutoff


class ToolError(Exception):
    """A tool refused to run (missing prerequisite, unknown locus, bad argument)."""


@dataclass
class ToolResult:
    tool: str
    headline: str
    text: str
    files: list[str] = field(default_factory=list)
    suggestions: list[str] = field(default_factory=list)

    def render(self) -> str:
        parts = [f"HEADLINE: {self.headline}", self.text.rstrip()]
        parts += [f"SUGGEST_FOLLOWUP: {s}" for s in self.suggestions]
        return "\n".join(p for p in parts if p)


@dataclass
class Sandbox:
    """Filesystem the tools work in. ``shared`` is the campaign-wide sandbox that
    later stages read (the paper's tasks shared one filesystem); ``artifacts`` is the
    current task's own record directory."""

    genomes: Path
    shared: Path
    artifacts: Path

    def contigs(self) -> dict[str, str]:
        return read_fasta(self.genomes)

    def table(self, name: str) -> list[dict[str, str]]:
        path = self.shared / name
        if not path.exists():
            raise ToolError(f"{name} does not exist yet; run the earlier stage first")
        lines = path.read_text().splitlines()
        cols = lines[0].split("\t")
        return [dict(zip(cols, ln.split("\t"), strict=True)) for ln in lines[1:] if ln]

    def publish(self, name: str, cols: list[str], rows: list[list]) -> str:
        body = "\t".join(cols) + "\n" + "".join(
            "\t".join(str(c) for c in r) + "\n" for r in rows)
        for base in (self.shared, self.artifacts):
            target = base / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(body)
        return name


# -- tool implementations ---------------------------------------------------------

def list_contigs(sb: Sandbox) -> ToolResult:
    """Stage 1: inventory the input genomes."""
    contigs = sb.contigs()
    rows = []
    for name, seq in contigs.items():
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        rows.append([name, len(seq), f"{gc:.3f}"])
    f = sb.publish("contigs.tsv", ["contig", "length_nt", "gc"], rows)
    total = sum(len(s) for s in contigs.values())
    return ToolResult("list_contigs", f"contigs={len(contigs)} total_nt={total}",
                      "\n".join(f"{r[0]}\t{r[1]} nt\tGC {r[2]}" for r in rows), [f])


def _translate(dna: str) -> str:
    return "".join(CODON_TABLE.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def _orfs(seq: str, min_aa: int) -> list[tuple[int, int, str]]:
    """(start, end, protein) of every ATG..stop ORF of >= min_aa residues, one strand."""
    found = []
    for frame in range(3):
        prot = _translate(seq[frame:])
        pos = 0
        for chunk in prot.split("*"):
            m = chunk.find("M")
            if m >= 0 and len(chunk) - m >= min_aa:
                s = frame + 3 * (pos + m)
                e = frame + 3 * (pos + len(chunk) + 1)  # include the stop codon
                found.append((s, e, chunk[m:]))
            pos += len(chunk) + 1
    return found


def find_rt_orfs(sb: Sandbox, min_aa: int = 100) -> ToolResult:
    """Stage 2: six-frame ORF search for the RT catalytic motif."""
    rows = []
    for contig, seq in sb.contigs().items():
        n = len(seq)
        for strand, s in (("+", seq), ("-", revcomp(seq))):
            for start, end, prot in _orfs(s, int(min_aa)):
                m = RT_MOTIF.search(prot)
                if not m:
                    continue
                fs, fe = (start, end) if strand == "+" else (n - end, n - start)
                rows.append([contig, strand, fs, fe, len(prot), m.group(0)])
    rows.sort(key=lambda r: (r[0], r[2]))
    rows = [[f"L{i:02d}", *r] for i, r in enumerate(rows, 1)]
    f = sb.publish("rt_candidates.tsv",
                   ["locus_id", "contig", "strand", "start", "end", "aa_len", "motif"], rows)
    return ToolResult("find_rt_orfs",
                      f"rt_orfs={len(rows)} contigs_with_rt={len({r[1] for r in rows})}",
                      "\n".join("\t".join(str(c) for c in r) for r in rows), [f])


def classify_rt(sb: Sandbox) -> ToolResult:
    """Stage 3: group RT candidates into classes, set aside fragments."""
    cands = sb.table("rt_candidates.tsv")
    label = {m: chr(ord("A") + i) for i, m in enumerate(sorted({c["motif"] for c in cands}))}
    rows, counts = [], {}
    for c in cands:
        frag = int(c["aa_len"]) < FRAGMENT_MAX_AA
        cls = "fragment" if frag else label[c["motif"]]
        counts[cls] = counts.get(cls, 0) + 1
        rows.append([c["locus_id"], cls, c["motif"], int(frag)])
    f = sb.publish("rt_classes.tsv", ["locus_id", "rt_class", "motif", "fragment"], rows)
    dist = ",".join(f"{k}:{v}" for k, v in sorted(counts.items()))
    return ToolResult("classify_rt",
                      f"classified={len(rows)} fragments={counts.get('fragment', 0)} "
                      f"classes={dist}",
                      "\n".join("\t".join(str(c) for c in r) for r in rows), [f])


def _upstream(sb: Sandbox, cand: dict[str, str]) -> str:
    seq = sb.contigs()[cand["contig"]]
    s, e = int(cand["start"]), int(cand["end"])
    if cand["strand"] == "+":
        return seq[max(0, s - MAX_UPSTREAM_SCAN * 2):s]
    return revcomp(seq[e:e + MAX_UPSTREAM_SCAN * 2])


def _rng(locus: str) -> random.Random:
    return random.Random(f"demo-{locus}")


@functools.lru_cache(maxsize=64)
def _scan(locus: str, upstream: str):
    # pure: the rng is seeded from the locus, so repeated calls (revisions, reruns)
    # may share one result
    return kmer_scan(locus, upstream, _rng(locus))


@functools.lru_cache(maxsize=64)
def _delimit(locus: str, upstream: str):
    return delimit_array(locus, upstream, _rng(locus))


def scan_arrays(sb: Sandbox) -> ToolResult:
    """Stage 4: paper array scan step 1 on the window upstream of every full-length RT."""
    cands = {c["locus_id"]: c for c in sb.table("rt_candidates.tsv")}
    classes = sb.table("rt_classes.tsv")
    rows, lines, suggestions = [], [], []
    for cl in classes:
        if cl["fragment"] == "1":
            continue
        lid = cl["locus_id"]
        call = _scan(lid, _upstream(sb, cands[lid]))
        rows.append([lid, call.status, call.R, call.seed or "-",
                     ",".join(map(str, call.copies)) or "-",
                     len(_upstream(sb, cands[lid]))])
        lines.append(f"{lid}\t{call.status}\tR={call.R}")
        if call.status == "array":
            suggestions.append(
                f"Deep dive locus {lid} (array, R={call.R}) using deep_dive locus={lid}")
    f = sb.publish("arrays.tsv", ["locus_id", "status", "R", "seed", "copies", "upstream_nt"],
                   rows)
    stat = [r[1] for r in rows]
    return ToolResult(
        "scan_arrays",
        f"assessed={len(rows)} arrays={stat.count('array')} "
        f"no_array={stat.count('no_array')} not_assessed={stat.count('not_assessed')}",
        "\n".join(lines), [f], suggestions)


def deep_dive(sb: Sandbox, locus: str = "all") -> ToolResult:
    """Stage 5: paper array delimitation step 2 plus spacing and RT-distance stats."""
    arrays = [r for r in sb.table("arrays.tsv") if r["status"] == "array"]
    done = {p.stem for p in (sb.shared / "deep_dives").glob("*.tsv")}
    if locus == "all":
        targets = [r for r in arrays if r["locus_id"] not in done]
    else:
        targets = [r for r in arrays if r["locus_id"] == locus]
        if not targets:
            raise ToolError(f"{locus} has no array call in arrays.tsv; nothing to dive into")
    cands = {c["locus_id"]: c for c in sb.table("rt_candidates.tsv")}
    cols = ["locus_id", "status", "copies", "block_len", "median_spacing_nt",
            "last_copy_to_rt_nt", "repeat_consensus"]
    lines, delimited, files = [], [], []
    for row in targets:
        lid = row["locus_id"]
        up = _upstream(sb, cands[lid])
        da = _delimit(lid, up)
        if da is None:
            rec = [lid, "not_delimited", "-", "-", "-", "-", "-"]
        else:
            window = up[-MAX_UPSTREAM_SCAN:]
            gaps = da.spacings or [b - a for a, b in zip(da.copy_starts, da.copy_starts[1:],
                                                         strict=False)]
            rec = [lid, "delimited", len(da.copy_starts), len(da.repeat),
                   int(statistics.median(gaps)), len(window) - da.copy_starts[-1], da.repeat]
            delimited.append(lid)
        files.append(sb.publish(f"deep_dives/{lid}.tsv", cols, [rec]))
        lines.append("\t".join(map(str, rec)))
    return ToolResult(
        "deep_dive",
        f"deep_dives={len(targets)} delimited={len(delimited)} "
        f"loci={','.join(r['locus_id'] for r in targets) or '-'}",
        "\n".join(lines), files)


TOOLS: dict[str, tuple[Callable[..., ToolResult], str]] = {
    "list_contigs": (list_contigs, "list_contigs  - inventory the input contigs"),
    "find_rt_orfs": (find_rt_orfs, "find_rt_orfs [min_aa=100]  - six-frame search for RT-motif ORFs"),
    "classify_rt": (classify_rt, "classify_rt  - class per RT candidate, set fragments aside"),
    "scan_arrays": (scan_arrays, "scan_arrays  - tandem-repeat scan of the window upstream of each RT"),
    "deep_dive": (deep_dive, "deep_dive [locus=<id>|all]  - delimit one array and measure it"),
}


def parse_call(arg_text: str, tool: str) -> dict[str, str]:
    """`k=v` tokens, or a single bare token bound to the tool's first parameter."""
    kwargs: dict[str, str] = {}
    for tok in arg_text.replace("`", " ").split():
        tok = tok.strip("[]")
        if not tok:
            continue
        if "=" in tok:
            k, v = tok.split("=", 1)
            kwargs[k.strip()] = v.strip()
        elif tool == "deep_dive":
            kwargs["locus"] = tok
        elif tool == "find_rt_orfs":
            kwargs["min_aa"] = tok
    return kwargs


def run_tool(sb: Sandbox, name: str, arg_text: str = "") -> ToolResult:
    if name not in TOOLS:
        raise ToolError(f"unknown tool {name!r}; available: {', '.join(TOOLS)}")
    fn = TOOLS[name][0]
    try:
        return fn(sb, **parse_call(arg_text, name))
    except (TypeError, ValueError) as exc:
        raise ToolError(f"bad arguments for {name}: {exc}") from exc
