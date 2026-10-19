# Folder validation: Python against original R gRodon

**PASS: 40 of 40 assemblies passed.** 0 numerical/settings/count mismatches; 0 errors; 0 not run. Inventory entries: 40.

## Inputs and comparison method

The supplied inventory contains 3 genomic DNA files and 37 protein-only files. Protein sequences cannot recover synonymous codons. Protein-only inputs require genomic DNA fetched for the same assembly accessions; the proteins are not back-translated. Each accession and both source/genome SHA256 identities are retained in the machine-readable report.

The input preparation audit recorded identical protein identifiers, sequences, and order for 37 of 37 protein inputs when compared with the current files for those exact NCBI assembly accessions. 32 compressed protein archives matched the official checksum; 5 differed. For 5 of the differing files, only descriptive annotation headers changed: the protein IDs, sequences, and order were identical. These are annotation text differences, not evidence of different assembly genomes. Accessions: `GCF_000010425.1`, `GCF_964063315.1`, `GCF_964063395.1`, `GCF_964063405.1`, `GCF_964084875.1`.

Both implementations consume the same predicted coding sequences and the same selected ribosomal gene IDs. Each row compares all twelve full-mode statistics, prepared gene counts, retained highly expressed gene counts, and prediction settings. A pass requires finite values and positive, ordered doubling-time confidence bounds. This checks agreement with the original implementation, not biological accuracy against measured growth rates.

Numerical acceptance: `abs(Python − R) <= max(1e-12, 1e-10 × max(abs(Python), abs(R)))`. The separate exact-equality count compares binary64 encodings after JSON parsing, including the sign of zero. Numerical agreement does not imply bitwise identity. R serializes its reported floating-point values to JSON, which can itself introduce final-digit rounding.

Every produced run artifact is checked against the SHA256 and byte-size manifest; missing, modified, or unrecorded artifacts fail the audit.

Inventory SHA256: `42542c4b894b3d4c789d3c5f0568d2354177de432b1d8f0e3917a338866987be`.

## Software and settings

- PYTHON: biopython 1.88; grodon2-py 0.2.0; numpy 2.5.4; pandas 3.0.6; pyrodigal 3.7.1; python 3.13.5; scipy 1.18.1
- R: Biostrings 2.80.2; R R version 4.6.1 (2026-06-24 ucrt); coRdon 1.30.0; gRodon 2.7.3; jsonlite 2.0.0; matrixStats 1.5.0
- BLASTP: 2.17.0

Observed pipeline settings (all distinct configurations). Recorded executable paths are historical provenance from the validation machine; use the portable reproduction commands below for a new checkout.

```json
{
  "blastp": "C:\\Users\\drgarza\\Documents\\GitHub\\gRodon2_py\\.reference\\managed-blast\\bin\\blastp.exe",
  "closed": true,
  "confidence": 0.95,
  "gene_mode": "single",
  "genetic_code": 11,
  "mask": true,
  "max_evalue": 1e-08,
  "min_identity": 70.0,
  "min_query_coverage": 0.8,
  "min_reference_coverage": 0.8,
  "min_ribosomal": 10,
  "mode": "full",
  "reference": "bacteria-swissprot-v1",
  "temperature": null,
  "threads": 1,
  "timeout": 1800,
  "training_set": "madin"
}
```

## Every assembly

Doubling times and their absolute difference are shown in minutes. `nHE` is the retained highly expressed gene count. The result column gives the run status, statistics within tolerance, and bitwise-equal statistics. Comparisons use unrounded values; only this display is rounded. Full confidence intervals and hour values are preserved in the JSON and CSV exports.

| Assembly | Organism | Input | CDS | nHE Python/R | Python min | R min | Absolute Δ min | Result (tolerance; exact) |
|---|---|---|---:|---:|---:|---:|---:|---|
| GCA_964062575.1 | Lactococcus lactis | genomic | 2832 | 41/41 | 7.27762771346 | 7.27762771346 | 3.16413562018e-14 | pass (12/12; 2/12) |
| GCA_964063365.1 | Rahnella aquatilis | genomic | 5497 | 44/44 | 26.3383228074 | 26.3383228074 | 2.6645352591e-14 | pass (12/12; 2/12) |
| GCA_964063435.1 | Leuconostoc mesenteroides | genomic | 2078 | 44/44 | 23.5830260634 | 23.5830260634 | 1.66533453694e-14 | pass (12/12; 2/12) |
| GCF_000010425.1 | Bifidobacterium adolescentis ATCC 15703 | protein | 1682 | 43/43 | 39.8056757981 | 39.8056757981 | 1.33226762955e-14 | pass (12/12; 3/12) |
| GCF_964060345.1 | Lentilactobacillus buchneri | protein | 2532 | 41/41 | 38.8834436045 | 38.8834436045 | 6.66133814775e-15 | pass (12/12; 4/12) |
| GCF_964062525.1 | Enterococcus gilvus | protein | 3533 | 42/42 | 17.7724668303 | 17.7724668303 | 3.66373598126e-14 | pass (12/12; 3/12) |
| GCF_964062605.1 | Lentilactobacillus kefiri | protein | 2515 | 41/41 | 40.8008961775 | 40.8008961775 | 6.66133814775e-15 | pass (12/12; 2/12) |
| GCF_964062635.1 | Pediococcus ethanolidurans | protein | 2194 | 41/41 | 118.995145178 | 118.995145178 | 3.46389583683e-13 | pass (12/12; 2/12) |
| GCF_964062735.1 | Lentilactobacillus parabuchneri | protein | 2496 | 41/41 | 41.5016748394 | 41.5016748394 | 6.66133814775e-14 | pass (12/12; 3/12) |
| GCF_964062745.1 | Levilactobacillus cerevisiae | protein | 3199 | 42/42 | 35.7442753488 | 35.7442753488 | 1.26565424807e-13 | pass (12/12; 2/12) |
| GCF_964062795.1 | Loigolactobacillus coryniformis | protein | 2927 | 42/42 | 70.4893896445 | 70.4893896445 | 5.3290705182e-14 | pass (12/12; 2/12) |
| GCF_964062975.1 | Lactococcus lactis | protein | 2294 | 42/42 | 13.1495027122 | 13.1495027122 | 2.6645352591e-14 | pass (12/12; 2/12) |
| GCF_964063155.1 | Pediococcus parvulus | protein | 2099 | 41/41 | 95.0539263342 | 95.0539263342 | 1.33226762955e-13 | pass (12/12; 2/12) |
| GCF_964063315.1 | Leuconostoc carnosum | protein | 1731 | 43/43 | 30.6229539717 | 30.6229539717 | 3.33066907388e-14 | pass (12/12; 2/12) |
| GCF_964063325.1 | Hafnia alvei | protein | 4626 | 44/44 | 42.0973441512 | 42.0973441512 | 5.99520433298e-14 | pass (12/12; 2/12) |
| GCF_964063335.1 | Hafnia alvei | protein | 4547 | 44/44 | 49.5306891246 | 49.5306891246 | 1.33226762955e-14 | pass (12/12; 4/12) |
| GCF_964063345.1 | Pseudomonas lini | protein | 5848 | 43/43 | 36.9859831962 | 36.9859831962 | 6.66133814775e-15 | pass (12/12; 2/12) |
| GCF_964063355.1 | Latilactobacillus curvatus | protein | 2177 | 42/42 | 20.9504332445 | 20.9504332445 | 0 | pass (12/12; 4/12) |
| GCF_964063375.1 | Bacillus pumilus | protein | 3844 | 43/43 | 20.7462873577 | 20.7462873577 | 9.99200722163e-15 | pass (12/12; 2/12) |
| GCF_964063385.1 | Latilactobacillus curvatus | protein | 2028 | 42/42 | 21.7493019544 | 21.7493019544 | 0 | pass (12/12; 3/12) |
| GCF_964063395.1 | Pediococcus parvulus | protein | 1836 | 41/41 | 103.411110282 | 103.411110282 | 2.93098878501e-13 | pass (12/12; 2/12) |
| GCF_964063405.1 | Leuconostoc mesenteroides | protein | 1907 | 44/44 | 23.5435871937 | 23.5435871937 | 6.66133814775e-15 | pass (12/12; 2/12) |
| GCF_964063415.1 | Latilactobacillus sakei | protein | 2072 | 42/42 | 20.3813966223 | 20.3813966223 | 6.66133814775e-15 | pass (12/12; 2/12) |
| GCF_964063425.1 | Lactiplantibacillus pentosus | protein | 3563 | 42/42 | 19.4250422867 | 19.4250422867 | 1.33226762955e-14 | pass (12/12; 2/12) |
| GCF_964063445.1 | Pediococcus parvulus | protein | 2195 | 41/41 | 97.9133936068 | 97.9133936068 | 3.99680288865e-14 | pass (12/12; 2/12) |
| GCF_964063755.1 | Latilactobacillus sakei | protein | 1919 | 42/42 | 20.3149258704 | 20.3149258704 | 2.6645352591e-14 | pass (12/12; 2/12) |
| GCF_964065095.1 | Lentilactobacillus buchneri | protein | 2425 | 41/41 | 38.4964400972 | 38.4964400972 | 1.13242748512e-13 | pass (12/12; 2/12) |
| GCF_964065115.1 | Levilactobacillus brevis | protein | 2430 | 42/42 | 48.6771169507 | 48.6771169507 | 2.6645352591e-14 | pass (12/12; 4/12) |
| GCF_964065165.1 | Levilactobacillus yonginensis | protein | 2658 | 42/42 | 45.540467792 | 45.540467792 | 6.66133814775e-15 | pass (12/12; 3/12) |
| GCF_964065205.1 | Lacticaseibacillus paracasei | protein | 3025 | 43/43 | 116.943589059 | 116.943589059 | 2.13162820728e-13 | pass (12/12; 3/12) |
| GCF_964065235.1 | Leuconostoc pseudomesenteroides | protein | 2148 | 43/43 | 28.0983047434 | 28.0983047434 | 7.32747196253e-14 | pass (12/12; 2/12) |
| GCF_964065245.1 | Lactiplantibacillus plantarum | protein | 3082 | 42/42 | 31.1757240693 | 31.1757240693 | 4.66293670343e-14 | pass (12/12; 2/12) |
| GCF_964065285.1 | Levilactobacillus brevis | protein | 2409 | 42/42 | 44.5705677827 | 44.5705677827 | 7.32747196253e-14 | pass (12/12; 2/12) |
| GCF_964065435.1 | Pediococcus parvulus | protein | 1911 | 41/41 | 96.9844479364 | 96.9844479364 | 5.3290705182e-14 | pass (12/12; 3/12) |
| GCF_964065445.1 | Lactiplantibacillus plantarum | protein | 3288 | 42/42 | 31.7113654588 | 31.7113654588 | 1.33226762955e-14 | pass (12/12; 2/12) |
| GCF_964065465.1 | Liquorilactobacillus satsumensis | protein | 2468 | 42/42 | 153.197196998 | 153.197196998 | 2.6645352591e-14 | pass (12/12; 3/12) |
| GCF_964068245.1 | Lentilactobacillus kefiri | protein | 2437 | 41/41 | 39.0097261716 | 39.0097261716 | 4.66293670343e-14 | pass (12/12; 3/12) |
| GCF_964069135.1 | Schleiferilactobacillus harbinensis | protein | 3032 | 41/41 | 98.2574787372 | 98.2574787372 | 3.46389583683e-13 | pass (12/12; 2/12) |
| GCF_964084875.1 | Leuconostoc carnosum | protein | 1724 | 43/43 | 30.2562335721 | 30.2562335721 | 5.99520433298e-14 | pass (12/12; 5/12) |
| GCF_964084915.1 | Lactiplantibacillus plantarum | protein | 2996 | 42/42 | 31.256114463 | 31.256114463 | 4.66293670343e-14 | pass (12/12; 2/12) |

## Agreement by statistic

| Statistic | Compared | Finite pairs | Within tolerance | Bitwise equal | Max absolute Δ | Max relative Δ |
|---|---:|---:|---:|---:|---:|---:|
| CUBHE | 40 | 40 | 40 | 3 | 5.55111512313e-16 | 6.80805433826e-16 |
| GC | 40 | 40 | 40 | 3 | 4.99600361081e-16 | 1.27387851258e-15 |
| GCdiv | 40 | 40 | 40 | 2 | 4.99600361081e-16 | 4.6340125563e-15 |
| ConsistencyHE | 40 | 40 | 40 | 1 | 5.55111512313e-16 | 1.05673775603e-15 |
| CUB | 40 | 40 | 40 | 6 | 6.66133814775e-16 | 1.18000501956e-15 |
| CPB | 40 | 40 | 40 | 2 | 4.4408920985e-16 | 1.45269499208e-15 |
| FilteredSequences | 40 | 40 | 40 | 40 | 0 | 0 |
| nHE | 40 | 40 | 40 | 40 | 0 | 0 |
| dCUB | 40 | 40 | 40 | 1 | 8.32667268469e-16 | 1.86765581428e-15 |
| d | 40 | 40 | 40 | 2 | 5.77315972805e-15 | 4.34775691305e-15 |
| LowerCI | 40 | 40 | 40 | 0 | 7.1054273576e-15 | 6.03448332629e-15 |
| UpperCI | 40 | 40 | 40 | 0 | 1.42108547152e-14 | 7.71762286816e-15 |

## Problems and warnings

No input, artifact, execution, or numerical comparison failures were found.

Warnings are retained even for passing comparisons:

- python: Filtered 113 of 1907 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063405.1`
- python: Filtered 114 of 1724 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964084875.1`
- python: Filtered 115 of 1731 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063315.1`
- python: Filtered 122 of 1836 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063395.1`
- python: Filtered 124 of 1919 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063755.1`
- python: Filtered 129 of 2468 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065465.1`
- python: Filtered 132 of 1911 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065435.1`
- python: Filtered 134 of 2148 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065235.1`
- python: Filtered 135 of 2099 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063155.1`
- python: Filtered 137 of 2028 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063385.1`
- python: Filtered 153 of 2194 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964062635.1`
- python: Filtered 154 of 2177 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063355.1`
- python: Filtered 156 of 2195 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063445.1`
- python: Filtered 156 of 2409 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065285.1`
- python: Filtered 157 of 2072 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063415.1`
- python: Filtered 168 of 2437 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964068245.1`
- python: Filtered 170 of 2294 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964062975.1`
- python: Filtered 172 of 2425 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065095.1`
- python: Filtered 174 of 2430 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065115.1`
- python: Filtered 175 of 2078 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCA_964063435.1`
- python: Filtered 180 of 3032 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964069135.1`
- python: Filtered 181 of 2496 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964062735.1`
- python: Filtered 185 of 2515 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964062605.1`
- python: Filtered 190 of 2658 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065165.1`
- python: Filtered 193 of 3533 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964062525.1`
- python: Filtered 197 of 2532 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964060345.1`
- python: Filtered 200 of 2927 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964062795.1`
- python: Filtered 201 of 2996 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964084915.1`
- python: Filtered 214 of 3025 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065205.1`
- python: Filtered 216 of 3082 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065245.1`
- python: Filtered 232 of 2832 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCA_964062575.1`
- python: Filtered 235 of 3199 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964062745.1`
- python: Filtered 271 of 4547 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063335.1`
- python: Filtered 282 of 3563 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063425.1`
- python: Filtered 284 of 4626 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063325.1`
- python: Filtered 286 of 3288 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964065445.1`
- python: Filtered 291 of 5848 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063345.1`
- python: Filtered 347 of 5497 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCA_964063365.1`
- python: Filtered 359 of 3844 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_964063375.1`
- python: Filtered 74 of 1682 genes: lengths must exceed 240 bp, be divisible by three, and contain only A/C/G/T — `GCF_000010425.1`
- r: Models were trained with a default genetic code of '1' for eukaryotes and '11' for prokaryotes. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi — `GCA_964062575.1`, `GCA_964063365.1`, `GCA_964063435.1`, `GCF_000010425.1`, `GCF_964060345.1`, `GCF_964062525.1`, `GCF_964062605.1`, `GCF_964062635.1`, `GCF_964062735.1`, `GCF_964062745.1`, `GCF_964062795.1`, `GCF_964062975.1`, `GCF_964063155.1`, `GCF_964063315.1`, `GCF_964063325.1`, `GCF_964063335.1`, `GCF_964063345.1`, `GCF_964063355.1`, `GCF_964063375.1`, `GCF_964063385.1`, `GCF_964063395.1`, `GCF_964063405.1`, `GCF_964063415.1`, `GCF_964063425.1`, `GCF_964063445.1`, `GCF_964063755.1`, `GCF_964065095.1`, `GCF_964065115.1`, `GCF_964065165.1`, `GCF_964065205.1`, `GCF_964065235.1`, `GCF_964065245.1`, `GCF_964065285.1`, `GCF_964065435.1`, `GCF_964065445.1`, `GCF_964065465.1`, `GCF_964068245.1`, `GCF_964069135.1`, `GCF_964084875.1`, `GCF_964084915.1`
- r: There were 1 genes with ambiguous bases, these genes have been ignored — `GCF_000010425.1`
- r: There were 113 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063405.1`
- r: There were 114 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964084875.1`
- r: There were 115 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063315.1`
- r: There were 122 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063395.1`
- r: There were 124 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063755.1`
- r: There were 129 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964065465.1`
- r: There were 132 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964065435.1`
- r: There were 134 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964065235.1`
- r: There were 135 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063155.1`
- r: There were 137 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063385.1`
- r: There were 153 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964062635.1`
- r: There were 154 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063355.1`
- r: There were 156 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063445.1`, `GCF_964065285.1`
- r: There were 157 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063415.1`
- r: There were 168 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964068245.1`
- r: There were 170 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964062975.1`
- r: There were 172 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964065095.1`
- r: There were 174 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964065115.1`
- r: There were 175 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCA_964063435.1`
- r: There were 180 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964069135.1`
- r: There were 181 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964062735.1`
- r: There were 185 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964062605.1`
- r: There were 190 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964065165.1`
- r: There were 193 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964062525.1`
- r: There were 197 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964060345.1`
- r: There were 200 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964062795.1`
- r: There were 201 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964084915.1`
- r: There were 214 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964065205.1`
- r: There were 216 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964065245.1`
- r: There were 232 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCA_964062575.1`
- r: There were 235 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964062745.1`
- r: There were 271 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063335.1`
- r: There were 282 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063425.1`
- r: There were 284 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063325.1`
- r: There were 286 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964065445.1`
- r: There were 291 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063345.1`
- r: There were 347 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCA_964063365.1`
- r: There were 359 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_964063375.1`
- r: There were 73 genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored — `GCF_000010425.1`

## Reproduce this test case

The archived report, input inventory, and compact reference fixture preserve their original contents and checksum relationships. Historical absolute paths are provenance, not required locations on another machine. Raw input genomes, the original protein folder, local tools, and complete run directories are not distributed. See the [benchmark archive notes](../validation/README.md).

The public reproduction below fetches the same forty genomic DNA files by their archived NCBI URLs and SHA256 checksums; it does not need the original source folder. The new inventory therefore contains forty DNA inputs, whereas the historical inventory records three DNA and thirty-seven protein inputs. The genomic sequence bytes used for prediction are the same. New reports are written separately from the archived evidence.

Run these PowerShell commands from the repository root with Python 3.11 or later. The genome extra pins Pyrodigal 3.7.1, and setup installs BLASTP 2.17.0 into `.tools`. For the full comparison, separately install original R gRodon 2.7.3, its dependencies, and jsonlite, then replace the example Rscript path with your installation. The observed versions above describe the archived runs.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev,genome]"
.\.venv\Scripts\python.exe -m grodon setup --tools-dir .tools/blast-2.17.0 --version 2.17.0
.\.venv\Scripts\python.exe scripts/fetch_benchmark_genomes.py --manifest validation/test_genomes_inputs_manifest.json --output-dir validation/test_genomes_inputs/public-genomes --fetch
.\.venv\Scripts\python.exe scripts/prepare_folder_genomes.py --input-dir validation/test_genomes_inputs/public-genomes --cache-dir validation/test_genomes_inputs/ncbi --manifest validation/reproduced_inputs.json --fetch
.\.venv\Scripts\python.exe scripts/validate_genome_folder.py --manifest validation/reproduced_inputs.json --runs-dir validation/test_genomes_runs --blastp .tools/blast-2.17.0/bin/blastp.exe --rscript "C:\path\to\R\bin\Rscript.exe" --workers 8
.\.venv\Scripts\python.exe scripts/report_folder_validation.py --manifest validation/reproduced_inputs.json --runs-dir validation/test_genomes_runs --output-json validation/reproduced_comparison.json --output-csv validation/reproduced_comparison.csv --output-md validation/reproduced_comparison.md
```

On Linux or macOS, use `python3 -m venv .venv`, `.venv/bin/python`, and `.tools/blast-2.17.0/bin/blastp`, with the appropriate Rscript path. For an offline repeat, omit `--fetch` from both genome fetching and preparation to verify the existing input cache. Add `--resume` to the runner to reuse completed runs only after their recorded inputs, settings, and artifact hashes have been checked. To regenerate every analysis, select a fresh `--runs-dir` for both the runner and reporter and omit `--resume`. After dependencies are installed, network access is confined to input preparation; prediction, R comparison, reporting, and cached-input verification run locally.

After a new complete run passes, its numeric R results can be exported separately. The second command below runs the existing archived model regression tests offline; it needs neither R nor the downloaded genomes:

```powershell
.\.venv\Scripts\python.exe scripts/export_folder_reference.py --report validation/reproduced_comparison.json --output validation/reproduced_growth_reference.json
.\.venv\Scripts\python.exe -m pytest tests/test_folder_growth_reference.py
```

That regression test compares fitted-model predictions with recorded original R results. It does not repeat gene prediction, BLAST, or codon-statistic computation; the full workflow above exercises those stages.

Repository artifacts: [full comparison JSON](../validation/test_genomes_comparison.json), [per-assembly CSV](../validation/test_genomes_comparison.csv), [input inventory](../validation/test_genomes_inputs_manifest.json), and [compact model regression fixture](../tests/data/folder_growth_reference.json). Input provenance follows the official [NCBI genome download documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/) and [NCBI Datasets REST API documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/).

Generated at 2026-10-19T10:53:21.421824+00:00. The JSON companion preserves every per-metric value, difference, exact-equality flag, setting, version, artifact hash, input identity check, and error.
