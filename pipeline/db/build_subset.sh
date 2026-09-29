#!/usr/bin/env bash
# pipeline/db/build_subset.sh -- build the filtered protein subset database.
#
# paper Methods "Metagenomic sequence database" p.28-29.
# paper: prodigal-gv 2.10.0, default parameters, for CDS calling.
# paper filters: exclude incomplete CDS, non-standard amino acids, proteins
# >8,000 residues, degenerate k-mer repeats, and >=50% low-complexity
# sequence (tantan).
# paper: DIAMOND linear-time clustering at 90% then 70% identity with >=80%
# coverage of the shorter sequence, then cascaded 50% identity with >=80%
# mutual coverage.
# paper: representative = member closest to the 80th percentile of cluster
# length.
#
# Usage: build_subset.sh [--dry-run] [GENOME_FNA] [OUTDIR]
# Defaults (env or positional):
#   GENOME_FNA  input genomes/contigs   (default data/db/genomes.fna)
#   OUTDIR      output dir              (default results/db)
#   ANNOTATED_FAA (env) source gene calls for annotated input assemblies —
#               retained verbatim per paper p.28 (prodigal-gv is for
#               unannotated assemblies only)
# NOT-IN-PAPER: paths and seqkit/tantan plumbing are ours.
set -euo pipefail

DRY_RUN=0
POS=()
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        -h|--help) grep '^#' "$0"; exit 0 ;;
        *) POS+=("$arg") ;;
    esac
done

GENOME_FNA="${POS[0]:-${GENOME_FNA:-data/db/genomes.fna}}"
OUT="${POS[1]:-${OUTDIR:-results/db}}"

# shellcheck source=../lib/run.sh
source "$(dirname "${BASH_SOURCE[0]}")/../lib/run.sh"

# Apply the paper's sequence filters; writes filtered FASTA.
# NOT-IN-PAPER: filter implementation is ours (stdlib python over prodigal-gv
# + tantan output); the paper states the criteria only.
# GAP: "degenerate k-mer repeats" -- k and the cutoff are not stated;
# placeholder: drop sequences whose distinct 3-mer count is <= 3 (effectively
# a single short repeat unit).
apply_filters() {
    python3 - "$1" "$2" "$3" <<'PY'
import sys

AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


masked = {h.split()[0]: seq for h, seq in read_fasta(sys.argv[2])}
# keys are FIRST TOKENS: tantan preserves the original header, and the lookup
# below keys on the first token — keying this dict on the FULL header made the
# get() always miss, so the low-complexity filter never fired (S1 finding B3)

seen_seqs = set()  # paper: "Proteins were deduplicated" (NOT-IN-PAPER: exact-
# sequence key, first occurrence kept; the paper does not state the method)

with open(sys.argv[3], "w") as out:
    for header, seq in read_fasta(sys.argv[1]):
        if "partial=00" not in header:  # paper: exclude incomplete CDS
            continue
        if not set(seq) <= AA20:  # paper: exclude non-standard amino acids
            continue
        if len(seq) > 8000:  # paper: exclude >8,000 residues
            continue
        kmers = {seq[i:i + 3] for i in range(len(seq) - 2)}
        if len(kmers) <= 3:  # paper: exclude degenerate k-mer repeats (GAP: k/cutoff)
            continue
        mask = masked.get(header.split()[0], "")
        low = sum(1 for c in mask if c.islower() or c == "X")
        if mask and low / len(mask) >= 0.5:  # paper: >=50% low-complexity (tantan)
            continue
        if seq in seen_seqs:  # paper: deduplicated
            continue
        seen_seqs.add(seq)
        out.write(f">{header}\n{seq}\n")
PY
}

# Representative = member closest to the 80th percentile of cluster length.
# Reads a rep<TAB>member cluster table + the member FASTA; writes rep ids.
# NOT-IN-PAPER: selection plumbing is ours.
select_reps() {
    python3 - "$1" "$2" "$3" <<'PY'
import sys

lengths = {}
name, n = None, 0
with open(sys.argv[2]) as fh:
    for line in fh:
        if line.startswith(">"):
            if name is not None:
                lengths[name] = n
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
if name is not None:
    lengths[name] = n

clusters = {}
with open(sys.argv[1]) as fh:
    for line in fh:
        rep, member = line.rstrip("\n").split("\t")[:2]
        clusters.setdefault(rep, []).append(member)

with open(sys.argv[3], "w") as out:
    for rep, members in clusters.items():
        ls = sorted(lengths.get(m, 0) for m in members)
        p80 = ls[min(len(ls) - 1, int(0.8 * (len(ls) - 1)))]
        best = min(members, key=lambda m: (abs(lengths.get(m, 0) - p80), m))
        out.write(best + "\n")
PY
}

if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "$OUT"
fi

# paper p.28 (verbatim): "Genes were predicted with prodigal-gv 2.10.0 (71,
# 77) under default parameters in unannotated assemblies of Logan, ENA, and
# NCBI, and the source gene calls were retained for annotated assemblies of
# NCBI and JGI." ANNOTATED_FAA (env) = existing source gene calls for the
# input assemblies: retain them verbatim; only unannotated input gets
# prodigal-gv 2.10.0 default parameters.
if [[ -n "${ANNOTATED_FAA:-}" ]]; then
    run_sh 'cp "$1" "$2"  # retain source gene calls (annotated assemblies)' \
        "$ANNOTATED_FAA" "$OUT/proteins_raw.faa"
else
    run prodigal-gv -i "$GENOME_FNA" -a "$OUT/proteins_raw.faa"
fi

# paper: >=50% low-complexity via tantan (tantan masks low-complexity as
# lowercase; -x X alternative noted)
run_sh 'tantan "$1" > "$2"' "$OUT/proteins_raw.faa" "$OUT/proteins_tantan.faa"

# paper filters: incomplete CDS, non-standard amino acids, >8,000 residues,
# degenerate k-mer repeats, >=50% low-complexity
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "python3 - $OUT/proteins_raw.faa $OUT/proteins_tantan.faa $OUT/proteins_filt.faa  # paper filters: complete CDS only; standard aa; <=8000 aa; non-degenerate k-mers; <50% low-complexity"
else
    apply_filters "$OUT/proteins_raw.faa" "$OUT/proteins_tantan.faa" "$OUT/proteins_filt.faa"
fi

# paper: DIAMOND linear-time clustering at 90% then 70% identity with >=80%
# coverage of the shorter sequence (cascaded on representatives), then 50%
# identity with >=80% mutual coverage.
# GAP: the paper does not give the DIAMOND version or exact coverage flag
# `diamond cluster` coverage flags are --member-cover / --mutual-cover (percents);
# --cov-mode is not a diamond cluster option (grok review 2026-09-24)
# coverage follow the mmseqs convention -- verify against the installed release.
run diamond makedb --in "$OUT/proteins_filt.faa" -d "$OUT/proteins_filt"
run diamond cluster -d "$OUT/proteins_filt.dmnd" -o "$OUT/clusters90.tsv" --approx-id 90 --member-cover 80
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "python3 - $OUT/clusters90.tsv $OUT/proteins_filt.faa $OUT/reps90.ids  # rep = member closest to 80th pct of length"
else
    select_reps "$OUT/clusters90.tsv" "$OUT/proteins_filt.faa" "$OUT/reps90.ids"
fi
run_sh 'seqkit grep -f "$1" "$2" > "$3"' "$OUT/reps90.ids" "$OUT/proteins_filt.faa" "$OUT/reps90.faa"
run diamond makedb --in "$OUT/reps90.faa" -d "$OUT/reps90"
run diamond cluster -d "$OUT/reps90.dmnd" -o "$OUT/clusters70.tsv" --approx-id 70 --member-cover 80
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "python3 - $OUT/clusters70.tsv $OUT/reps90.faa $OUT/reps70.ids  # rep = member closest to 80th pct of length"
else
    select_reps "$OUT/clusters70.tsv" "$OUT/reps90.faa" "$OUT/reps70.ids"
fi
run_sh 'seqkit grep -f "$1" "$2" > "$3"' "$OUT/reps70.ids" "$OUT/reps90.faa" "$OUT/reps70.faa"
run diamond makedb --in "$OUT/reps70.faa" -d "$OUT/reps70"
run diamond cluster -d "$OUT/reps70.dmnd" -o "$OUT/clusters50.tsv" --approx-id 50 --mutual-cover 80
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "python3 - $OUT/clusters50.tsv $OUT/reps70.faa $OUT/subset_reps.ids  # final subset representatives (80th pct of length)"
else
    select_reps "$OUT/clusters50.tsv" "$OUT/reps70.faa" "$OUT/subset_reps.ids"
fi
run_sh 'seqkit grep -f "$1" "$2" > "$3"' "$OUT/subset_reps.ids" "$OUT/proteins_filt.faa" "$OUT/subset.faa"
