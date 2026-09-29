# Real-genome array scan — re-run receipt (2026-09-29)

Verification-ladder step 3 (REPRODUCTION.md §2): `scripts/scan_genome.py` on real phage
genomes must reproduce the MarsHill 5-copy shuffle-controlled array call.

## Inputs
Fetched from ENA on 2026-09-29 (`https://www.ebi.ac.uk/ena/browser/api/fasta/<acc>`):

| accession | phage | length (nt) |
|---|---|---|
| MW248466.1 | Staphylococcus phage MarsHill | 266,637 |
| MW218148.1 | Staphylococcus phage vB_StaM_SA1 | 260,727 |
| MW349129.1 | Staphylococcus phage Madawaska | 265,446 |

## Command
```
python3 scripts/scan_genome.py MW248466.1.fna MW218148.1.fna MW349129.1.fna
```
Run on `main` @ 0f3f150, Python 3.13.5, fixed seed 20260923 (in the script), 142 s.
Output: `scan-main-0f3f150.txt` (the script prints the top 8 arrays per genome).

## Result
- **Anchor reproduced.** MarsHill `pos 213,994-214,801  R=5 copies  span=0.81 kb
  seed=ATATGAATACGTAT  shuffle_max=0`. The same seed gives the 5-copy call in SA1
  (9,597-10,371) and Madawaska (213,195-214,001).
- **Array counts differ from the 2026-09-24 receipt**
  (`task-runs/20260924-devin-build/real-genome-scan.txt`):

  | genome | 2026-09-24 receipt | main @ 0f3f150 |
  |---|---|---|
  | MW248466.1 | 36 | 50 |
  | MW218148.1 | 16 | 26 |
  | MW349129.1 | 36 | 57 |

## Why the counts differ (verified)
The seed is fixed, so this is code drift, not randomness. Running the script at the
commit that produced the old receipt (93a7b47) on MW248466.1 gives 36 arrays again
(`scan-MarsHill-at-93a7b47.txt`), the old number exactly. The scan and `arrays.py`
were changed by later review rounds (window-local suppression and coordinates, skipped-copy
chaining, tolerant windows), and the genome-wide call list grew.

## Not established
- Whether the additional calls are true arrays or false positives. The only ground truth
  here is the MarsHill anchor; the paper reports arrays at 28/95 ART loci, not genome-wide
  phage scans, so the 50/26/57 totals are not comparable to any paper number.
- This is the standalone genome-wide driver, not the per-locus scan the paper ran
  upstream of each RT (needs RT coordinates and the ART family from the full pipeline).
