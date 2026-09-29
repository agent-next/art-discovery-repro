#!/usr/bin/env bash
# pipeline/rnaseq/sa1_infection.sh -- SA1 infection time-course RNA-seq.
#
# paper Methods p.36-37: BioProject PRJNA836150 -- SA1 infection of
# Staphylococcus lentus, 12 libraries over 3 time points.
# paper: fastp 1.3.6, infection runs state ONLY "minimum length 30 nt"
# (Methods p.36). The old adapter-off/poly-G/min-12 flags here belonged to the
# SMALL-RNA plasmid libraries (a different dataset, p.36-37) -- S1 finder B,
# verified against the PDF text.
# paper: Bowtie2 --very-sensitive -X 1000 --no-unal against the SA1 genome
# (MW218148.1) plus the S. lentus chromosome (NZ_CP059679.1).
# paper: keep MAPQ >= 10; TPM over 259 features.
#
# Usage: sa1_infection.sh [--dry-run] [ACCESSIONS]
#   ACCESSIONS  file with one SRA run accession per line (the 12 libraries)
#               (default data/rnaseq/PRJNA836150.accessions)
# Env vars (defaults are ours; NOT-IN-PAPER path plumbing):
#   REF_SA1     SA1 genome FASTA        (paper: MW218148.1)
#   REF_HOST    S. lentus chromosome    (paper: NZ_CP059679.1)
#   FEATURES    feature annotation GTF  (paper: 259 features)
#   OUTDIR      output dir              (default results/rnaseq)
#   THREADS     worker threads for the tools that take them (default 1; NOT-IN-PAPER)
#   FASTQ_SOURCE  sra (default: prefetch + fasterq-dump, ~10 GB of plain FASTQ per
#               library) or ena (the run's .fastq.gz from ENA, ~1.4 GB); NOT-IN-PAPER
#   KEEP_INTERMEDIATES  1 (default) keeps FASTQ/trimmed/SAM/unsorted BAM; 0 deletes them
#               after each library's sorted BAM exists (12 libraries are ~17 GB
#               compressed, ~10x that as FASTQ + SAM; NOT-IN-PAPER disk plumbing)
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

ACCESSIONS="${POS[0]:-${ACCESSIONS:-data/rnaseq/PRJNA836150.accessions}}"
REF_SA1="${REF_SA1:-data/rnaseq/MW218148.1.fna}"
REF_HOST="${REF_HOST:-data/rnaseq/NZ_CP059679.1.fna}"
FEATURES="${FEATURES:-data/rnaseq/features_259.gtf}"
OUT="${OUTDIR:-results/rnaseq}"
THREADS="${THREADS:-1}"
FASTQ_SOURCE="${FASTQ_SOURCE:-sra}"

# SECURITY (S3b 2026-09-24): argument form -- never eval. Data-derived values
# (accessions from a file) cross a trust boundary here.
run() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf '%s\n' "$*"
    else
        "$@"
    fi
}

# SRA-style accession gate: the only file-derived input allowed near a command.
acc_ok() { [[ "$1" =~ ^[A-Z]{3}[A-Za-z0-9_-]{0,20}$ ]]; }  # shell-metachar-free

# TPM from featureCounts output: TPM = (count/len) / sum(count/len) * 1e6.
# paper: TPM over 259 features. NOT-IN-PAPER: counter/normalizer plumbing.
tpm_calc() {
    python3 - "$1" "$2" <<'PY'
import sys

header = None
rows = []
with open(sys.argv[1]) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "Geneid":
            header = f
            continue
        rows.append(f)
rpk = [[float(c) / float(f[5]) for c in f[6:]] for f in rows]
scale = [s / 1e6 for s in (sum(col) for col in zip(*rpk))]
with open(sys.argv[2], "w") as out:
    out.write("Geneid\t" + "\t".join(header[6:]) + "\n")
    for f, col in zip(rows, rpk):
        tpm = [c / s if s else 0.0 for c, s in zip(col, scale)]
        out.write(f[0] + "\t" + "\t".join(f"{v:.6g}" for v in tpm) + "\n")
PY
}

if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "$OUT/fastq" "$OUT/trim" "$OUT/bam"
fi

# paper: combined reference of SA1 genome MW218148.1 + S. lentus NZ_CP059679.1
# FP-filter follow-up (2026-09-24): positional bash -c -- no interpolation
# into the command string at all.
run bash -c 'cat "$1" "$2" > "$3"' bash "$REF_SA1" "$REF_HOST" "$OUT/sa1_plus_host.fna"
run bowtie2-build --threads "$THREADS" "$OUT/sa1_plus_host.fna" "$OUT/bt2_sa1_host"

# NOT-IN-PAPER: accession-file plumbing. In --dry-run without the file we still
# print the per-library commands against 12 placeholder run ids (paper: 12
# libraries), so the dry run is self-contained.
acc_stream() {
    if [[ -f "$ACCESSIONS" ]]; then
        cat "$ACCESSIONS"
    elif [[ "$DRY_RUN" -eq 1 ]]; then
        printf 'SRR_PLACEHOLDER_%02d\n' $(seq 1 12)
    else
        echo "missing accessions file: $ACCESSIONS" >&2
        return 1
    fi
}

while IFS= read -r acc; do
    [[ -z "$acc" || "$acc" == \#* ]] && continue
    if ! acc_ok "$acc"; then
        printf 'skipping non-accession line: %s\n' "$acc" >&2
        continue
    fi
    # paper: BioProject PRJNA836150, 12 paired-end libraries
    # NOT-IN-PAPER: fetch/prefetch plumbing
    if [[ "$FASTQ_SOURCE" == ena ]]; then
        run bash -c 'set -o pipefail
            for n in 1 2; do
                url=$(curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp&format=tsv" |
                    tail -n +2 | cut -f2 | tr ";" "\n" | grep "_${n}.fastq.gz")
                curl -fsSL -o "$2/$1_${n}.fastq.gz" "https://${url}"
            done' bash "$acc" "$OUT/fastq"
        r1="$OUT/fastq/${acc}_1.fastq.gz" r2="$OUT/fastq/${acc}_2.fastq.gz"
    else
        run prefetch -O "$OUT/fastq" "$acc"
        run fasterq-dump --split-files -e "$THREADS" -O "$OUT/fastq" "$acc"
        r1="$OUT/fastq/${acc}_1.fastq" r2="$OUT/fastq/${acc}_2.fastq"
    fi
    # paper: "fastp 1.3.6 (minimum length 30 nt)". Adapter/poly-G behavior is
    # unstated for the infection runs -- NOT-IN-PAPER: fastp defaults apply.
    run fastp --length_required 30 --thread "$THREADS" \
        -i "$r1" -I "$r2" \
        -o "$OUT/trim/${acc}_1.fq.gz" -O "$OUT/trim/${acc}_2.fq.gz"
    # paper: Bowtie2 --very-sensitive -X 1000 --no-unal vs SA1 + host
    # paper: "Properly paired alignments with MAPQ of at least 10 and a template
    # of at most 1,500 nt were retained as fragments." Flags verified against the
    # htslib samtools-view manual: -f/--require-flags (0x2 = proper pair),
    # -e/--expr with the documented `tlen` variable.
    # NOT-IN-PAPER plumbing: the aligner streams into the filter instead of
    # writing an intermediate SAM (~13 GB per library here); same tools, same flags.
    run bash -c 'set -o pipefail
        bowtie2 --very-sensitive -X 1000 --no-unal -x "$1" -p "$2" -1 "$3" -2 "$4" |
        samtools view -b -q 10 -f 2 -e "tlen <= 1500 && tlen >= -1500" -@ "$2" -o "$5" -' \
        bash "$OUT/bt2_sa1_host" "$THREADS" "$OUT/trim/${acc}_1.fq.gz" \
        "$OUT/trim/${acc}_2.fq.gz" "$OUT/bam/${acc}.bam"
    run samtools sort -o "$OUT/bam/${acc}.sorted.bam" -@ "$THREADS" "$OUT/bam/${acc}.bam"
    if [[ "${KEEP_INTERMEDIATES:-1}" == 0 ]]; then
        run rm -rf "$OUT/fastq/${acc}" "$r1" "$r2" \
            "$OUT/trim/${acc}_1.fq.gz" "$OUT/trim/${acc}_2.fq.gz" \
            "$OUT/bam/${acc}.bam"
    fi
done < <(acc_stream)

# paper: "Each fragment was counted once on the strand of read 2 (the sense
# read)"; features = 258 annotated SA1 CDS + the array RNA (259 total).
# Subread Users Guide (verified): -p counts fragments; "For paired-end reads,
# strand of the first read is taken as the strand of the whole fragment", so
# -s 2 (reversely stranded) = the strand of READ 2 on proper pairs.
# -t CDS matches CDS-style annotation (default 'exon' would silently give 0
# counts). GAP: the paper's midpoint-assignment rule needs a custom counter;
# featureCounts assigns by overlap (NOT-IN-PAPER tool choice).
run featureCounts -p -s 2 -t CDS -T "$THREADS" -a "$FEATURES" -o "$OUT/counts.tsv" "$OUT"/bam/*.sorted.bam
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "python3 - $OUT/counts.tsv $OUT/tpm.tsv  # TPM = (count/len)/sum(count/len)*1e6"
else
    tpm_calc "$OUT/counts.tsv" "$OUT/tpm.tsv"
fi
