set -uo pipefail
export MAMBA_ROOT_PREFIX=~/tools/mamba
eval "$(~/tools/mm/bin/micromamba shell hook -s bash)"; micromamba activate art
cd ~/runs/real
export THREADS=8 KEEP_INTERMEDIATES=0 OUTDIR=results/rnaseq REF_SA1=data/rnaseq/MW218148.1.fna REF_HOST=data/rnaseq/NZ_CP059679.1.fna FEATURES=data/rnaseq/features.gtf
for acc in $(echo SRR19152326); do
  [ -f results/rnaseq/bam/$acc.sorted.bam ] && continue
  while [ "$(df --output=avail -BG / | tail -1 | tr -dc 0-9)" -lt 12 ]; do echo "disk low, waiting"; sleep 60; done
  echo "$acc" > data/rnaseq/one.accession
  bash ~/runs/art-harness-20260929/pipeline/rnaseq/sa1_infection.sh data/rnaseq/one.accession || echo "FAILED $acc"
  echo "LIB_DONE $acc $(date +%T)"
done
echo RNASEQ_DONE
