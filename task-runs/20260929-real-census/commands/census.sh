set -euo pipefail
REPO=${REPO:?repo checkout}; DATA=${DATA:?data root}; cd "$DATA"; mkdir -p results/census
micromamba activate art   # env with HMMER 3.4, MMseqs2, seqkit
zcat data/seq/viral.1.protein.faa.gz > data/seq/seqdb.faa
echo "proteins: $(grep -c ">" data/seq/seqdb.faa)"
time bash $REPO/pipeline/census/01_hmmsearch.sh data/hmm/rt_pfam6.hmm data/seq/seqdb.faa results/census/hmmsearch.domtblout
grep -vc "^#" results/census/hmmsearch.domtblout || true
time python3 $REPO/pipeline/census/02_filter.py results/census/hmmsearch.domtblout --seqs data/seq/seqdb.faa --out results/census/retained.ids
wc -l results/census/retained.ids
seqkit grep -f results/census/retained.ids data/seq/seqdb.faa > results/census/filtered.faa
grep -c ">" results/census/filtered.faa
time bash $REPO/pipeline/census/03_cluster.sh results/census/filtered.faa results/census/clu results/census/tmp
grep -c ">" results/census/clu_rep_seq.fasta
echo CENSUS_DONE
