#!/usr/bin/env bash
# pipeline/art_family/phylogeny.sh -- RT-domain phylogeny of the ART family.
#
# paper Methods "Human-directed analyses" p.33: 774-sequence set (98 ART+NART,
# 617 retron, 59 other RTs).
# paper: RT domain taken as the Pfam RVT_1 envelope +/- 40 residues.
# paper: MAFFT L-INS-i alignment; trimal -gappyout -> 254 columns.
# paper: IQ-TREE 3.1.2, ModelFinder (BIC) -> Q.pfam+F+R6; 1,000 ultrafast
# bootstrap + 1,000 SH-aLRT replicates; midpoint rooted.
#
# Usage: phylogeny.sh [--dry-run] [INPUT_FAA] [OUTDIR]
# Defaults (env or positional):
#   INPUT_FAA  the 774-sequence set      (default data/art_family/rt_774.faa)
#   RVT1_HMM   Pfam RVT_1 profile HMM    (default data/hmm/RVT_1.hmm)
#   OUTDIR     output dir                (default results/art_family/phylogeny)
# NOT-IN-PAPER: paths and seqkit/gotree plumbing are ours.
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

INPUT_FAA="${POS[0]:-${INPUT_FAA:-data/art_family/rt_774.faa}}"
OUT="${POS[1]:-${OUTDIR:-results/art_family/phylogeny}}"
RVT1_HMM="${RVT1_HMM:-data/hmm/RVT_1.hmm}"

# shellcheck source=../lib/run.sh
source "$(dirname "${BASH_SOURCE[0]}")/../lib/run.sh"

if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "$OUT"
fi

# paper: RT domain = Pfam RVT_1 envelope +/- 40 residues
# NOT-IN-PAPER: envelope extraction plumbing (domtblout env coords 20/21, bed
# is 0-based half-open).
run hmmsearch --noali --domtblout "$OUT/rvt1.domtbl" -o /dev/null "$RVT1_HMM" "$INPUT_FAA"
run_sh "$(cat <<'EOS'
awk '!/^#/ {s=$20-40; if (s<1) s=1; print $1"\t"s-1"\t"$21+40}' "$1" | sort -u > "$2"
EOS
)" "$OUT/rvt1.domtbl" "$OUT/rvt1_env.bed"
run_sh 'seqkit subseq --bed "$1" "$2" > "$3"' "$OUT/rvt1_env.bed" "$INPUT_FAA" "$OUT/rt_domain.faa"

# paper: MAFFT L-INS-i (--localpair --maxiterate 1000)
run_sh 'mafft --localpair --maxiterate 1000 "$1" > "$2"' "$OUT/rt_domain.faa" "$OUT/rt_domain.aln.faa"

# paper: trimal -gappyout -> 254 columns
run trimal -gappyout -in "$OUT/rt_domain.aln.faa" -out "$OUT/rt_domain.trim.faa"

# paper: IQ-TREE 3.1.2; ModelFinder with BIC (-m MFP) selected Q.pfam+F+R6;
# 1,000 ultrafast bootstrap (-B 1000) + 1,000 SH-aLRT (-alrt 1000)
run iqtree3 -s "$OUT/rt_domain.trim.faa" -m MFP -B 1000 -alrt 1000 -T AUTO --prefix "$OUT/rt_tree"

# paper: midpoint rooted
# NOT-IN-PAPER: gotree for the rooting step (paper states midpoint only)
run gotree reroot midpoint -i "$OUT/rt_tree.treefile" -o "$OUT/rt_tree.midpoint.nwk"
