#!/usr/bin/env bash
# pipeline/art_family/family_definition.sh -- ART family delineation.
#
# paper chain (Methods "Human-directed analyses" p.31-32):
#   * BLASTP seeded from the MarsHill RT (GenBank QQM14740.1)
#   * profile HMM from 12 GenBank phage RTs (10 Staphylococcus + 2 Listeria
#     phage LPJP1) via mafft --auto + hmmbuild
#   * hmmsearch (HMMER 3.4) E <= 1e-5 -> 9,100 proteins; keep 350-900 aa -> n=4,379
#   * pool with the 12 phage RTs + 2,019 myRT reference RTs
#   * diamond blastp --very-sensitive E <= 1e-5; edges at bitscore >= 140;
#     connected component containing the MarsHill RT = 117 proteins (no
#     reference RT)
#   * expand clusters, add members >= 400 aa (706 further) -> 823 full-length RTs
#   * mmseqs easy-cluster --min-seq-id 0.9 -c 0.8 -> 230 clusters; the member on
#     the longest contig represents each cluster
#   * FAMSA alignment + FastTree LG
#   * ART clade = smallest clade containing every locus with a type I partner
#     and >=500 nt of non-coding upstream; 95 representatives ART_01..ART_95
#   * geNomad 1.12.0, database release 1.9, end-to-end default settings
#
# Usage: family_definition.sh [--dry-run]
# Inputs via env vars (defaults are ours; NOT-IN-PAPER path plumbing):
#   SEED_FAA      MarsHill RT FASTA (paper: GenBank QQM14740.1)
#   SEED_ID       accession inside SEED_FAA            (default QQM14740.1)
#   PHAGE_RT_FAA  12 GenBank phage RTs (paper: 10 Staphylococcus + 2 LPJP1)
#   MYRT_REF_FAA  2,019 myRT reference RTs (paper)
#   SEARCHDB_FAA  protein search database (metagenomic database representatives)
#   GENBANK_BLASTDB  BLAST protein database of GenBank (makeblastdb output prefix)
#   CONTIG_LEN_TSV   protein_id<TAB>contig length for every clustered member
#   CLUST_TSV     precomputed cluster table (rep<TAB>member) for the expansion
#   GENOMES_FNA   genomes/contigs for geNomad
#   GENOMAD_DB    geNomad database dir (paper: release 1.9)
#   OUTDIR        output dir (default results/art_family)
set -euo pipefail

DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        -h|--help) grep '^#' "$0"; exit 0 ;;
        *) echo "unknown arg: $arg" >&2; exit 2 ;;
    esac
done

SEED_FAA="${SEED_FAA:-data/art_family/QQM14740.1.faa}"
SEED_ID="${SEED_ID:-QQM14740.1}"
PHAGE_RT_FAA="${PHAGE_RT_FAA:-data/art_family/phage_rt_12.faa}"
MYRT_REF_FAA="${MYRT_REF_FAA:-data/art_family/myrt_ref_2019.faa}"
SEARCHDB_FAA="${SEARCHDB_FAA:-data/seq/seqdb.faa}"
GENBANK_BLASTDB="${GENBANK_BLASTDB:-data/blastdb/genbank_protein}"
CONTIG_LEN_TSV="${CONTIG_LEN_TSV:-data/art_family/contig_len.tsv}"
CLUST_TSV="${CLUST_TSV:-results/census/clu_cluster.tsv}"
GENOMES_FNA="${GENOMES_FNA:-data/art_family/genomes.fna}"
GENOMAD_DB="${GENOMAD_DB:-data/genomad_db_1.9}"
OUT="${OUTDIR:-results/art_family}"

# shellcheck source=../lib/run.sh
source "$(dirname "${BASH_SOURCE[0]}")/../lib/run.sh"

# BFS over DIAMOND self-hits: connected component containing SEED_ID.
# paper: edges at bitscore >= 140; component = 117 proteins, no reference RT.
# NOT-IN-PAPER: component extraction implementation is ours (stdlib python).
component_ids() {
    python3 - "$1" "$2" "$3" <<'PY'
import sys

edges = {}
with open(sys.argv[1]) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12 or float(f[11]) < 140.0:
            continue
        edges.setdefault(f[0], set()).add(f[1])
        edges.setdefault(f[1], set()).add(f[0])
seen = {sys.argv[2]}
stack = [sys.argv[2]]
while stack:
    for member in edges.get(stack.pop(), ()):
        if member not in seen:
            seen.add(member)
            stack.append(member)
with open(sys.argv[3], "w") as out:
    out.write("\n".join(sorted(seen)) + "\n")
PY
}

if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "$OUT"
fi

# --- seed BLASTP ----------------------------------------------------------
# paper: homologs of the MarsHill RT (GenBank QQM14740.1) were collected FROM
# GENBANK with BLASTP and, separately, from the metagenomic database with the
# profile HMM below. The paper's pool (4,379 database hits + 12 phage RTs + 2,019
# myRT references) does not contain the BLASTP hits: this branch supplies the
# GenBank phage RTs (the 12 profile sequences, and the "cultured phage" class).
# GAP: BLASTP parameters are not stated; -evalue 1e-5 mirrors the paper's
# hmmsearch threshold below, outfmt 6 plumbing is ours.
run blastp -query "$SEED_FAA" -db "$GENBANK_BLASTDB" -evalue 1e-5 -outfmt 6 -out "$OUT/seed_blastp.tsv"

# --- phage-RT profile HMM ---------------------------------------------------
# paper: profile HMM from 12 GenBank phage RTs (10 Staphylococcus + 2 Listeria
# phage LPJP1) via mafft --auto + hmmbuild
run_sh 'mafft --auto "$1" > "$2"' "$PHAGE_RT_FAA" "$OUT/phage_rt.aln.faa"
run_sh 'hmmbuild "$1" "$2" > "$3"' "$OUT/phage_rt.hmm" "$OUT/phage_rt.aln.faa" "$OUT/hmmbuild.log"

# --- hmmsearch --------------------------------------------------------------
# paper: hmmsearch (HMMER 3.4) E <= 1e-5 -> 9,100 proteins
run hmmsearch -E 1e-5 --tblout "$OUT/phage_rt.tbl" --domtblout "$OUT/phage_rt.domtbl" -o "$OUT/phage_rt.hmmsearch.txt" "$OUT/phage_rt.hmm" "$SEARCHDB_FAA"
# NOT-IN-PAPER: hit-id extraction plumbing (tblout column 1)
run_sh "$(cat <<'EOS'
grep -v '^#' "$1" | awk '{print $1}' | sort -u > "$2"
EOS
)" "$OUT/phage_rt.tbl" "$OUT/hits.ids"

# --- length filter ------------------------------------------------------------
# paper: keep 350-900 aa -> n=4,379
# NOT-IN-PAPER: seqkit for id/length plumbing (paper states thresholds only)
run_sh 'seqkit grep -f "$1" "$2" | seqkit seq -m 350 -M 900 > "$3"' "$OUT/hits.ids" "$SEARCHDB_FAA" "$OUT/hits_350_900.faa"

# --- sequence pool --------------------------------------------------------------
# paper: pool with the 12 phage RTs + 2,019 myRT reference RTs
run_sh 'cat "$1" "$2" "$3" > "$4"' "$OUT/hits_350_900.faa" "$PHAGE_RT_FAA" "$MYRT_REF_FAA" "$OUT/pool.faa"

# --- all-vs-all DIAMOND ---------------------------------------------------------
# paper: diamond blastp --very-sensitive E <= 1e-5; edges at bitscore >= 140
run diamond makedb --in "$OUT/pool.faa" -d "$OUT/pool"
run diamond blastp --very-sensitive -e 1e-5 -q "$OUT/pool.faa" -d "$OUT/pool" -o "$OUT/pool_self.tsv" --outfmt 6

# --- MarsHill connected component -------------------------------------------------
# paper: connected component containing the MarsHill RT = 117 proteins
# (no reference RT)
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "python3 - $OUT/pool_self.tsv $SEED_ID $OUT/component117.ids  # BFS component, edges bitscore >= 140"
else
    component_ids "$OUT/pool_self.tsv" "$SEED_ID" "$OUT/component117.ids"
fi
run_sh 'seqkit grep -f "$1" "$2" > "$3"' "$OUT/component117.ids" "$OUT/pool.faa" "$OUT/component117.faa"

# --- cluster expansion -------------------------------------------------------------
# paper: expand clusters, add members >= 400 aa (706 further) -> 823 full-length RTs
# NOT-IN-PAPER: expansion plumbing (rep<TAB>member table + seqkit length filter)
run_sh "$(cat <<'EOS'
awk 'NR==FNR {keep[$1]=1; next} keep[$1] {print $2}' "$1" "$2" | sort -u > "$3"
EOS
)" "$OUT/component117.ids" "$CLUST_TSV" "$OUT/expanded.ids"
run_sh 'seqkit grep -f "$1" "$2" | seqkit seq -m 400 > "$3"' "$OUT/expanded.ids" "$SEARCHDB_FAA" "$OUT/expanded_400.faa"
run_sh 'cat "$1" "$2" | seqkit rmdup -s > "$3"' "$OUT/component117.faa" "$OUT/expanded_400.faa" "$OUT/art_823.faa"

# --- 90% clustering --------------------------------------------------------------------
# paper: mmseqs easy-cluster --min-seq-id 0.9 -c 0.8 -> 230 clusters; "the member
# on the longest contig was taken as the representative of each cluster"
# (Methods p.31), replacing MMseqs' own representative.
run mmseqs easy-cluster "$OUT/art_823.faa" "$OUT/art_823_clu" "$OUT/tmp_mmseqs" --min-seq-id 0.9 -c 0.8
run python3 "$(dirname "${BASH_SOURCE[0]}")/longest_contig_reps.py" "$OUT/art_823_clu_cluster.tsv" "$CONTIG_LEN_TSV" "$OUT/art_reps.ids"
run_sh 'seqkit grep -f "$1" "$2" > "$3"' "$OUT/art_reps.ids" "$OUT/art_823.faa" "$OUT/art_reps.faa"

# --- representative alignment + tree -----------------------------------------------------
# paper: FAMSA alignment + FastTree LG (over the 230 cluster representatives)
run famsa "$OUT/art_reps.faa" "$OUT/art_reps.famsa.aln"
run_sh 'FastTree -lg "$1" > "$2"' "$OUT/art_reps.famsa.aln" "$OUT/art_reps.nwk"

# --- ART clade -----------------------------------------------------------------------------
# paper: ART clade = smallest clade containing every locus with a type I
# partner and >=500 nt of non-coding upstream; 95 representatives ART_01..ART_95
# GAP: the paper states the clade criterion but no delineation tool; applied to
# art_reps.nwk with locus context (phylogenetic curation step, not scripted).
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "# ART clade: smallest clade in $OUT/art_reps.nwk containing every locus with a type I partner + >=500 nt non-coding upstream -> ART_01..ART_95 (GAP: delineation tool not scripted)"
fi

# --- geNomad --------------------------------------------------------------------------------
# paper: geNomad 1.12.0, database release 1.9, end-to-end default settings
run genomad end-to-end "$GENOMES_FNA" "$OUT/genomad" "$GENOMAD_DB"
