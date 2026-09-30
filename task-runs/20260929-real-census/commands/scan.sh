export MAMBA_ROOT_PREFIX=~/tools/mamba
eval "$(~/tools/mm/bin/micromamba shell hook -s bash)"; micromamba activate art
cd ~/runs/art-harness-20260929
time python3 scripts/scan_rt_loci.py ~/runs/real/results/census/retained2.ids ~/runs/real/data/seq/viral.1.genomic.gbff.gz -o ~/runs/real/results/rt_loci.tsv -j 6
echo SCAN_DONE
