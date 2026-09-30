# Real-data runs on agent-dev-01 (2026-09-29)

Host: `agent-dev-01` (16 vCPU, 62 GB, shared 309 GB disk, other lanes use 91-100% of it).
Tools from micromamba env `art`: HMMER 3.4, MMseqs2, Bowtie2, SAMtools, subread
(featureCounts), fastp 1.3.7, sra-tools 3.4.1. Repo at `main` after PRs #24-#27.

Every small output is committed next to this file; only the multi-GB inputs (RefSeq
protein/GenBank flat files, FASTQ, BAMs) stay on the box and are re-downloadable.

| path | content |
| --- | --- |
| `commands/` | the exact wrapper scripts run on the box (`census.sh`, `rnaseq.sh`, `rnaseq326.sh`, `scan.sh`) |
| `census/` | `hmmsearch.domtblout`, `retained.ids` (271), `retained2.ids` (271, after the PR #24 fix), `filtered.faa`, `clu_rep_seq.fasta` (117), `clu_cluster.tsv` |
| `controls/` | control sequences (UniProt) and their domtblout / retained lists |
| `rnaseq/` | `features.gtf` (260 rows), `PRJNA836150.accessions`, `tpm_all.tsv` (fragments and TPM per feature for the 3 libraries run) |
| `rt_loci/rt_loci.tsv` | per-locus scan result (254 rows) |

Pfam profiles used (public, download from InterPro): PF00078 RVT_1, PF07727 RVT_2,
PF13456 RVT_3, PF13655 RVT_N, PF08388 GIIM, PF01348 Intron_maturas2.

## 1. Census subset (steps 01-03) on RefSeq viral proteins

Input: NCBI RefSeq `viral.1.protein.faa.gz`, 722,738 proteins. Six public Pfam RT
profiles (RVT_1, RVT_2, RVT_3, RVT_N, GIIM, Intron_maturas2), not the paper's 45 myRT
profiles (GAP, see REPRODUCTION.md).

```
hmmsearch --noali -Z 8 --domZ 8 -E 0.01 --domE 0.01 rt_pfam6.hmm seqdb.faa
```

Result: 772 hit proteins -> 271 pass the step-02 filter -> 117 clusters (step 03).

### Controls (paper reference RTs fetched from UniProt)

| control | retained before fix | retained after fix |
| --- | --- | --- |
| Ec86 (P23070) | yes | yes |
| BPP-1 Brt (Q775D8) | yes | yes |
| LtrA (P0A3U0) | **no** (best-domain coverage 0.68 < 0.75) | yes |
| AbiK (Q48614) | no | no (needs the myRT class profiles, not available) |
| RT-Cas1 fusion | not fetched | not fetched |

Retained controls in the 17-sequence control set: 12 -> 14. LtrA's RT core is split
over two domain hits; coverage is now the union of a profile's domains (PR #24, test
uses the real LtrA domtblout lines and failed before the fix).

## 2. RNA-seq reanalysis PRJNA836150 (SA1 + S. lentus)

Bowtie2 `--very-sensitive -X 1000 --no-unal`; samtools `-q 10 -f 2` with
`tlen <= 1500`; `featureCounts -p -s 2 -t CDS` on a GTF of the 258 SA1 CDS plus the array
RNA at 9,448-10,646 (built by `pipeline/rnaseq/genbank_to_gtf.py`, one feature per
strand, so 260 rows). Share = TPM of the array RNA over all phage features.

| library | condition | fragments assigned | array RNA, sense (+) | array RNA, antisense | array share of phage TPM |
| --- | --- | --- | --- | --- | --- |
| SRR19152324 | B+P_55 | 19,769,643 | 298,497 | 100 | 1.08% |
| SRR19152325 | B+P_15 | 8,648,798 | 704,175 | 374 | 4.52% |
| SRR19152327 | B+P_15 | 11,879,973 | 850,141 | 390 | 3.79% |

Paper: array-derived RNAs = 8% of intracellular phage RNA at 15 min.

Status: **not matched.** Two 15-min replicates give 4.5% and 3.8%, about half the
paper's value. The strand is right (sense:antisense about 2000:1) and the 15 min > 55
min ordering is right. Possible causes, none tested: the paper assigns a fragment by its
midpoint (featureCounts assigns by overlap), the paper's array-RNA feature may be
delimited differently, or the paper's 8% is a pooled or different-replicate figure.
SRR19152326 (third 15-min replicate) and the other 8 libraries were not run.

## 3. RT-locus scan with gene annotation (GAP-9)

`scripts/scan_rt_loci.py` (PR #27) on the 271 retained viral RTs against RefSeq viral
GenBank records: 254 loci scanned, 17 ids not in the flat file, 147 no_array, 103
not_assessed (upstream window too short), 4 array calls. Only one delimited to copies
(12 x GACTTGTT, an entomopoxvirus RT, RT-adjacent). This subset is eukaryotic-virus RTs,
so no phage RT-arrays are expected; the run shows the annotation rules execute on real
records, not that new arrays exist.

## 4. Disk incident

The shared disk reached 0 bytes free twice during the RNA-seq run (once from a 12.8 GB
SAM, once from other lanes while `fasterq-dump` was writing). Mitigations merged: bowtie2
streamed into samtools (no SAM, PR #25), `KEEP_INTERMEDIATES=0`, plain FASTQ deleted right
after fastp (PR #26). SRR19152326 failed with "storage exhausted" and was abandoned; a
disk-safe run of the remaining libraries needs a dedicated volume.

## Cleanup

Scratch removed from the box: FASTQ, trimmed reads, partial BAMs, tmux sessions. Kept for
now: `~/runs/real/results/rnaseq/bam/*.sorted.bam` (3 libraries) and the census results.
