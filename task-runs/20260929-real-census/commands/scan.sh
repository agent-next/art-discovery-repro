REPO=${REPO:?repo checkout}; DATA=${DATA:?data root}; cd "$REPO"
micromamba activate art
time python3 scripts/scan_rt_loci.py "$DATA/results/census/retained2.ids" "$DATA/data/seq/viral.1.genomic.gbff.gz" -o "$DATA/results/rt_loci.tsv" -j 6
echo SCAN_DONE
