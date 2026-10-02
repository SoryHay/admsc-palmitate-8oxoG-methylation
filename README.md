# admsc-palmitate-8oxoG-methylation

Analysis code for **"Lipotoxic Palmitate Hyperpolarizes Mitochondria, Raises 8-oxo-dG Damage in Depleted mtDNA, and Converts Nuclear
5mC to 5hmC Genome-Wide in Human Adipose-Derived Mesenchymal Stem Cells"** (Gospodinova A., Velichkov A., Dimitrov N., Todorova K.,
Hayrabedyan S.; International Journal of Molecular Sciences, submitted 2026).

Processed data: Zenodo https://doi.org/10.5281/zenodo.23111729 (region-level 5mC/5hmC and 8-oxo-dG profiles, figure source data).
Raw sequencing data are not public: they contain the genetic information of the cell donor (Regulation (EU) 2016/679, Art. 9) and
are available from the corresponding authors on reasonable request under a data use agreement.

## What is here
The scripts are archived **as they were run**; the working comments of the analysis were removed and each file carries a one-line description, the code itself is unchanged. Paths of the analysis server are kept, except machine-specific prefixes, which are
replaced by placeholders (`<DATA_DISK>`, `<NAS_SHARE>`, `<SESSION_TMP>`, `~`). `PROVENANCE.tsv` lists the original location of every file.
GPU steps ran in Docker containers pinned to one card.

| module | content | figures / tables |
|---|---|---|
| `01_basecalling_and_conversion` | Guppy + Rerio basecalling, modkit tag update, alignment to T2T-CHM13v2.0, modkit pileup, CpG/CpH split from the reference; dorado guide basecall and esox 8-oxo-dG calling; region annotation | all |
| `02_global_and_local_methylation` | pooled 5mC per region class and gene panel, per-read methylation, repeat classes, metagene, IL-1β axis, PoreMeth2 segmentation | Fig. 9, 11, 12a, S6, S9 |
| `03_5mC_to_5hmC` | 10 % read subset, dorado 5mCG_5hmCG basecall, pileup, Δ5mC/Δ5hmC, read-quality checks, ribosomal DNA | Fig. 13a–b, Tables S3–S5 |
| `04_8oxodG_and_5mC_around_lesions` | esox call filtering and read-space rates per compartment, position in read, sequence context, molecule methylation tertiles, 5mC around 8-oxo-dG anchors | Fig. 8, 10, 13c, S4, S5, S8 |
| `05_figures` | drawing scripts of Figures 8–13 and S5, S6, S9 | |
| `06_image_analysis` | `.ims` → OME-TIFF, Nellie network segmentation, Cellpose cell masks, per-cell and paired-cell analysis of the 24 h time-lapse | Fig. 5, 6, S1–S3 |
| `07_public_profiles` | builds the public Zenodo tables from the analysis tables | Zenodo record |

## Pipeline in pseudocode
```
# 5mC (Guppy/Rerio)                                        module 01, 02
for arm in NT, PA:
    guppy_basecaller 5.0.16  --config rerio res_dna_r941_min_modbases-all-context_v001  fast5_pass/  -> unaligned BAM (Mm/Ml)
    modkit update-tags --mode implicit (legacy Mm/Ml -> MM/ML);  samtools fastq -T MM,ML | minimap2 -y (T2T-CHM13v2.0)  -> aln_T2T.bam
    modkit pileup aln_T2T.bam --ref T2T --filter-threshold C:0.70 --filter-threshold A:0.70   -> bedMethyl
    keep rows whose reference base is C (+) / G (-); CpG if ref[p+1]=G (+) or ref[p-1]=C (-)  -> ctx_thr070.cpg.tsv.gz
    region % 5mC = sum(Nmod) / sum(Nvalid) over the region   (pooled fraction; regions with >= 20 valid calls in both arms)
    modkit extract calls --cpg  -> per-read calls -> per-read beta = modified / valid CpG calls of the read

# 5mC + 5hmC (dorado combined model)                        module 03
seeded 10 % of pass read IDs -> fast5_subset -> pod5 -> dorado 0.9.6 dna_r9.4.1_e8_sup@v3.3 + 5mCG_5hmCG@v0
-> minimap2 -> modkit pileup --cpg --modified-bases 5mC 5hmC --filter-threshold C:0.7
-> per arm: 5mC % and 5hmC % of valid CpG calls; Δ = PA - NT; repeated for Q >= 12, Q >= 14, quality deciles, rDNA interval sets

# 8-oxo-dG (esox)                                          module 01, 04
dorado 0.9.6 SUP guide basecall (FASTQ, bare @read_id) -> esox resquiggle -> bonito -> remora modcall (Docker)
keep guanines whose 5-mer is on the 110-mer allowlist; call = score >= per-k-mer threshold (esox calibration, never retuned)
membership of a read in a compartment from its primary alignment (chrM, gene panel, rDNA intervals)
rate = calls / callable guanines x 10^6, read space, no end trimming; Poisson / exact binomial intervals

# 5mC around 8-oxo-dG on the same read                     module 04
anchors: guanines with an esox score (>= 0.70 -> 8-oxo-dG anchor; otherwise ordinary guanine)
for every CpG 11-1000 bp from an anchor on the same read (±10 bp masked):  pair = (distance, Rerio 5mC prediction score)
classes: 8-oxo-dG anchors · ordinary G on lesion-carrying reads (canon_all - canon_free, d4) · ordinary G on lesion-free reads
mean score per 50 bp bin and class
```

## Not redistributed
- PoreMeth2 helper scripts `ModkitResorter.sh` and `ParseModkit.pl` (authors of PoreMeth2: Mattei et al., Genome Res 2025) — obtain
  them from the PoreMeth2 authors.
- `thresholds.json` (110 5-mer allowlist and per-k-mer thresholds) — from the esox publication (Pagès-Gallego et al., Nat Commun
  2025, doi:10.1038/s41467-025-60391-3).
- Reference: T2T-CHM13v2.0 (`01_basecalling_and_conversion/fetch_t2t.sh`).

## Notes
- `layer1_cpg_aggregate.py` also writes a per-region Fisher q-value column; it assumes independent CpG observations, which does not
  hold for nanopore reads, and is not used in the article. Only the pooled fractions and ratios are.
- One library per arm (n = 1): the code reports pooled estimates, not tests between conditions.

## Software
Guppy 5.0.16 · Rerio all-context v001 · dorado 0.9.6 · esox (Docker image) · modkit 0.6.4 · minimap2 2.29 · samtools · bedtools ·
pysam · numpy · pandas · scipy · matplotlib · R + PoreMeth2 · methylartist 1.5.2 · Python 3.11, Nellie 1.0.4, napari 0.8.0,
Cellpose 4.2.1.1, laptrack 0.17.1.

## Licence
MIT (see `LICENSE`).
