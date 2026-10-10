set -uo pipefail
REPO=${REPO:?repo checkout}; DATA=${DATA:?data root}; cd "$DATA"
micromamba activate art   # env with Bowtie2, samtools, subread, fastp, sra-tools
export THREADS=8 KEEP_INTERMEDIATES=0 OUTDIR=results/rnaseq REF_SA1=data/rnaseq/MW218148.1.fna REF_HOST=data/rnaseq/NZ_CP059679.1.fna FEATURES=data/rnaseq/features.gtf
for acc in $(cat data/rnaseq/rest11.accessions); do
  [ -f results/rnaseq/bam/$acc.sorted.bam ] && continue
  while [ "$(df --output=avail -BG / | tail -1 | tr -dc 0-9)" -lt 17 ]; do echo "disk low, waiting"; sleep 60; done
  echo "$acc" > data/rnaseq/one.accession
  bash "$REPO/pipeline/rnaseq/sa1_infection.sh" data/rnaseq/one.accession || echo "FAILED $acc"
  echo "LIB_DONE $acc $(date +%T)"
done
echo RNASEQ_DONE
