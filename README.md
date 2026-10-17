# Extensive H3K4me1 defines broad domains safeguarded against DNA methylation in mouse, cow and human oocytes

Analysis and figure code for:

> Fosslie M\*, Hashim A\*, *et al.* Extensive H3K4me1 defines broad domains safeguarded against DNA methylation in mouse, cow and human oocytes. (2026)
> \* equal contribution

**Data:** GEO accession `GSE332907`. Public datasets: Supp. Table 1 of the manuscript.

The scripts produce the data panels of the figures. Labels, colours, dashed threshold lines and panel
layout were finalized in Adobe Illustrator, so the raw script output differs from the published figures
in appearance, not in data. Every script also writes the plotted values as CSV/TXT.

---

## Figure index

| Figure | Script | Runs on | Input |
|--------|--------|---------|-------|
| Genome browser tracks (Fig. 1a, 3d, 4a, 4c; Supp. 1a, 1f, 2a, 2f, 4a, 4b) | `01_preprocessing/genomic_track_UCSC.sh` | server | BAM |
| Fig. 1b, 2b (heatmaps at broad domains) | `02_heatmaps/heatmap_deeptools.sh` | server | BigWig, peaks |
| Fig. 1c (genome coverage) | `03_genome_coverage/genome_coverage.sh` → `plot_genome_coverage.R` | server → RStudio | peaks |
| Fig. 1d (genomic elements) | `04_genomic_elements/genomic_element_distribution.R` | RStudio | peaks |
| Fig. 1e, 1f (H3K4me1 vs H3K27ac) | `05_mixture_model/Fig1e_1f_K4me1_vs_K27ac.R` | RStudio | `data/Fig1/` |
| Fig. 2f (clustering H3K36me3 / H3K4me1 / RNA, k = 15) | `06_kmeans/Fig2f_kmeans_K36me3_K4me1_RNA.R` | RStudio | `data/matrices/` |
| Fig. 3a–c (MII H3K4me3/H3K27me3 groups, 5 kb) | `05_mixture_model/Fig3abc_mouse_grouped_heatmap.R` | RStudio | `data/matrices/` |
| Supp. Fig. 3b (clustering, k = 5) | `06_kmeans/SuppFig3b_kmeans_MII_K4me3_K27me3.R` | RStudio | `data/matrices/` |
| Fig. 4b (bovine groups) | `05_mixture_model/Fig4b_bovine_grouped_heatmap.R` | RStudio | `data/Fig4/` |
| Fig. 4d (human groups) | `05_mixture_model/Fig4d_human_grouped_heatmap.R` | RStudio | `data/Fig4/` |
| Fig. 4e, Supp. Fig. 4c (human subclustering, 400 bp) | `06_kmeans/Fig4e_human_subclustering_400bp.R` | RStudio | `data/Fig4/` + Fig. 4d output |

The mouse bin matrices (`data/matrices/`) are made on the server with `07_bin_matrices/` (section 7) and
shared between figures: Fig. 2f and Supp. Fig. 3b read the same 10 kb matrix, Fig. 3a–c the 5 kb matrices.
Run Fig. 4d before Fig. 4e (4e uses the H3K4me3+/H3K27me3− bins written by 4d).

## Repository structure

```
.
├── README.md
├── environment.yml                    # conda environment with all software versions
├── H3K4me1-oocyte-analysis.Rproj      # open in RStudio: all paths relative to the repository root
├── 00_shared/
│   └── normalization.R                # 95th-percentile, top-5% and ln normalization; matrix reader
├── 01_preprocessing/
│   └── genomic_track_UCSC.sh          # BAM -> RPKM BigWig for genome browser tracks
├── 02_heatmaps/
│   └── heatmap_deeptools.sh           # size-sorted heatmaps at broad domains (deepTools)
├── 03_genome_coverage/
│   ├── genome_coverage.sh             # % genome covered by broad peaks (bedtools)
│   ├── plot_genome_coverage.R         # bar plot (Fig. 1c)
│   └── genome_coverage_labels.tsv     # sample -> plot label and order
├── 04_genomic_elements/
│   └── genomic_element_distribution.R # peaks over genomic elements (Fig. 1d)
├── 05_mixture_model/
│   ├── mixture_model.R                # shared mixture-model functions
│   ├── Fig1e_1f_K4me1_vs_K27ac.R      # Fig. 1e/1f
│   ├── Fig3abc_mouse_grouped_heatmap.R # Fig. 3a-c
│   ├── Fig4b_bovine_grouped_heatmap.R # Fig. 4b
│   └── Fig4d_human_grouped_heatmap.R  # Fig. 4d
├── 06_kmeans/
│   ├── kmeans_clustering.R            # shared k-means / heatmap functions
│   ├── Fig2f_kmeans_K36me3_K4me1_RNA.R      # Fig. 2f
│   ├── SuppFig3b_kmeans_MII_K4me3_K27me3.R  # Supp. Fig. 3b
│   └── Fig4e_human_subclustering_400bp.R    # Fig. 4e and Supp. Fig. 4c
├── 07_bin_matrices/
│   ├── make_bins.sh                   # genome-wide 5 kb / 10 kb bins (bedtools makewindows)
│   ├── build_chip_matrix.sh           # BAMs -> read-count matrix (bedtools multicov)
│   ├── build_methylation_matrix.sh    # BigWigs -> methylation matrix (multiBigwigSummary), 0-1
│   ├── bismark_cov_to_bigwig.sh       # Bismark coverage files -> methylation BigWig
│   ├── mm10_noMY.chrom.sizes          # mm10 chr1-19, chrX
│   └── samples/                       # sample sheets (matrix column label -> file)
├── peaks/                             # broad peak files (*.broadPeak) for Fig. 1c/1d (GEO; not tracked)
├── data/                              # input matrices, see data/README.md (not tracked)
│   ├── Fig1/                          #   Fig. 1e/1f
│   ├── Fig4/                          #   Fig. 4b, 4d, 4e, Supp. Fig. 4c
│   └── matrices/                      #   mouse 5 kb / 10 kb matrices (Fig. 2f, 3a-c, Supp. Fig. 3b)
├── results/                           # server output used by R scripts, e.g. genome_coverage/ (not tracked)
├── figures/                           # figure output, one folder per figure (created when run; not tracked)
└── logs/                              # SLURM logs (not tracked)
```

## Installation

```bash
git clone https://github.com/Adnanhashim/H3K4me1-oocyte-analysis.git
cd H3K4me1-oocyte-analysis
conda env create -f environment.yml
conda activate h3k4me1-oocyte
```

Server scripts (`*.sh`) run on a SLURM cluster (`sbatch`; set `--account` in the `#SBATCH` header first)
or on any Linux machine (`bash`; the `#SBATCH` lines are then ignored). Environment modules are loaded
when available (names as on the Saga HPC cluster); otherwise the conda environment is used. Submit from
the repository root after `mkdir -p logs`.

Figure scripts (`*.R`) run locally in RStudio: open `H3K4me1-oocyte-analysis.Rproj`, open the script and
click **Source**. Input paths are set at the top of each script. R packages:

```r
install.packages(c("ggplot2", "svglite", "hexbin", "gridExtra", "pheatmap", "RColorBrewer",
                   "mixtools", "fBasics", "data.table"))
BiocManager::install(c("GenomicRanges", "ChIPpeakAnno", "TxDb.Mmusculus.UCSC.mm10.knownGene"))
```

Each figure script writes its outputs to `figures/<figure>/`, including `sessionInfo.txt`.

## Software versions

Server analyses were run on the Saga HPC cluster (Sigma2 / NRIS, Norway):

| Tool        | Version | Used for                                   |
|-------------|---------|--------------------------------------------|
| Trim Galore | 0.6.7   | Adapter/quality trimming                   |
| FastQC      | 0.11.9  | Read quality control                       |
| MultiQC     | 1.11    | QC report aggregation                      |
| Bowtie2     | 2.4.4   | Read alignment                             |
| SAMtools    | 1.13    | BAM processing, indexing                   |
| HTSlib      | 1.12    | BAM/BGZF support                           |
| BamTools    | 2.5.2   | BAM processing                             |
| BEDTools    | 2.30.0  | Interval operations, bins, read counts     |
| Subread     | 2.0.3   | Read counting (featureCounts)              |
| deepTools   | 3.5.0   | BigWigs, heatmaps, methylation per bin     |
| R           | ≥ 4.0   | Mixture models, k-means, plotting          |

## Analysis steps

### 1. Genome browser tracks

`01_preprocessing/genomic_track_UCSC.sh` converts deduplicated BAM files to RPKM-normalized BigWig tracks
(10 bp bins; chrX and chrM excluded from normalization; duplicates ignored; reads extended to fragment length).

| Genome  | Effective genome size |
|---------|-----------------------|
| mm10    | 2,652,783,500         |
| hg38    | 2,913,022,398         |
| bosTau8 | 2,640,168,518         |

```bash
sbatch 01_preprocessing/genomic_track_UCSC.sh -g mm10 -o bigwig/ bam/*_dedup.bam
```

### 2. Heatmaps at broad domains (Fig. 1b, 2b)

`02_heatmaps/heatmap_deeptools.sh` computes signal (RPKM BigWig) in a ±15 kb window centred on each broad
domain (10 bp bins, blacklist excluded) and plots a white–black heatmap with domains sorted by size.
The values matrix (`*_values.tab`) and the region order (`*_sorted_regions.bed`) are exported as text.

```bash
sbatch 02_heatmaps/heatmap_deeptools.sh \
    -r peaks/P10_GO_H3K4me1.broadPeak -s bigwig/P10_GO_H3K4me1_rep1_dedup_RPKM.bw \
    -n P10_GO_H3K4me1_rep1 -l mm10-blacklist.v2.bed -z 22 -o heatmaps/
```

Options: `-z` colour maximum (default 22), `-w` window in bp (default 15000), `-L` regions label,
`-t` threads (`-h` for all options).

### 3. Genome coverage of broad peaks (Fig. 1c)

`03_genome_coverage/genome_coverage.sh` runs `bedtools coverage` of each peak file against the genome and
reports genome coverage as the mean per-chromosome covered fraction × 100 (`genome_coverage.txt`).
`plot_genome_coverage.R` plots the bar chart; labels and order are set in `genome_coverage_labels.tsv`
(tab-separated).

```bash
bash 03_genome_coverage/genome_coverage.sh -g mm10_sizes.bed -o results/genome_coverage/ peaks/*.broadPeak
```

Then copy `genome_coverage.txt` to `results/genome_coverage/` and Source `plot_genome_coverage.R`.
Output: `figures/Fig1/Fig1c_genome_coverage.svg` and the plotted values (CSV).

### 4. Genomic element distribution (Fig. 1d)

`04_genomic_elements/genomic_element_distribution.R` annotates every `*.broadPeak` file in `peaks/` with
`ChIPpeakAnno::genomicElementDistribution` (TxDb.Mmusculus.UCSC.mm10.knownGene; promoter −2000/+500 bp of
TSS; downstream 0–2000 bp of gene end): one EPS per sample and the values as CSV.

### 5. Mixture modelling: grouped heatmaps and Fig. 1e/1f

Bins are classified as positive or negative per dataset with two- or three-component Gaussian mixture
models (`05_mixture_model/mixture_model.R`, `mixtools::normalmixEM`, fixed seed, quantile starting
values). Each figure script sets its own parameters:

| Figure | Bins | Normalization | k | Threshold |
|--------|------|---------------|---|-----------|
| 1e/1f | mm10 10 kb, TSS-distal (> 2 kb from refTSS v3.3) | ÷ 95th percentile per dataset | 2 | mean of component means, fitted on ln values |
| 3a–c | mm10 5 kb (526,765 bins) | top 5% per dataset | – | fixed: MII H3K4me3 0.488, MII H3K27me3 0.042 (mean of 2 replicates)¹ |
| 4b | bosTau8 10 kb | top 5% per dataset | 3 | mean of the two lowest component means |
| 4d | hg38 10 kb | top 5% per dataset | 2 | mean of component means (values of the original fit) |

¹ Thresholds that reproduce the group sizes of the published Fig. 3a (34.9 / 13.0 / 19.3 / 32.8 % of
bins). A two-component fit on these 5 kb data is not stable (it can split off the low-signal tail), so
the fit is only reported (`thresholds.csv`). Set `fixed_thresholds <- NULL` in the script to use the fit.

Grouped heatmaps: bins are grouped as H3K27me3+ / overlap / H3K4me3+ / negative (random seeded order
within a group), shown as 0–1 white–black heatmaps with the group sizes (n) and, for Fig. 4b/4d, the
binary calls. Outputs: heatmap(s), group sizes, thresholds, mixture-model fits (PDF) and per-bin values.

Fig. 1e/1f: 2D histograms of H3K4me1 vs H3K27ac in TSS-distal 10 kb bins with Pearson r and the
mixture-model thresholds (dashed lines), and the percentage of H3K4me1-positive bins
that are H3K27ac-negative per cell type (PDF/SVG + CSV).

### 6. k-means clustering

All clustered heatmaps use `06_kmeans/kmeans_clustering.R`: read counts per bin (`bedtools multicov`)
normalized per dataset to the top 5% signal, bins without signal removed, `kmeans` (nstart = 25,
seed 123), clusters ordered by a fixed rule, 0–1 white–black heatmap (`pheatmap`) with cluster sizes (n).

| Figure | Bins | Clustered on | k | Cluster order |
|--------|------|--------------|---|---------------|
| 2f | mm10 10 kb | P10 GO, P14 GO, FGO H3K36me3 (2 replicates) and RNA-seq; P10 GO and 4w NSN H3K4me1 (2 replicates). DNA methylation (P10, P15, FGO) displayed only | 15 | H3K36me3-marked clusters by RNA, then H3K36me3-low clusters by H3K4me1 (as in the manuscript) |
| Supp. 3b | mm10 10 kb | P10 GO H3K4me1, MII H3K4me3, MII H3K27me3 (2 replicates each) | 5 | by size |
| 4e, Supp. 4c | 25 × 400 bp per H3K4me3+/H3K27me3− 10 kb bin (hg38, from Fig. 4d) | MI/II H3K4me1, 4-cell H3K4me3 (rep1, rep2), MII DNA methylation | 10 | by profile pattern (as in the manuscript) |

Outputs: heatmap (TIFF), bin → cluster table, cluster sizes and mean signal per dataset (CSV).

### 7. Bin matrices (server)

Bin-based figures read one read-count matrix and one DNA methylation matrix per genome and bin size.
Column names are the labels in the sample sheets (`07_bin_matrices/samples/`), which the figure scripts
use. Bins: mm10 chr1–19 and chrX (`mm10_noMY.chrom.sizes`): 5 kb → 526,765 bins (Fig. 3),
10 kb → 263,389 bins (Fig. 2f, Supp. Fig. 3b).

```bash
mkdir -p logs data/bins

# bins
bash 07_bin_matrices/make_bins.sh -g 07_bin_matrices/mm10_noMY.chrom.sizes -w 10000 -o data/bins/mm10_10kb_bins.bed
bash 07_bin_matrices/make_bins.sh -g 07_bin_matrices/mm10_noMY.chrom.sizes -w 5000  -o data/bins/mm10_5kb_bins.bed

# methylation data available only as Bismark coverage files (P10 GO, 2 replicates merged)
sbatch 07_bin_matrices/bismark_cov_to_bigwig.sh -g 07_bin_matrices/mm10_noMY.chrom.sizes \
       -o 5mC/P10_GO_GSM2690945_46_merged.bw 5mC/GSM2690945_*.cov.txt 5mC/GSM2690946_*.cov.txt

# read counts (ChIP-seq, picoChIP-seq, RNA-seq; BAMs counted in parallel)
sbatch 07_bin_matrices/build_chip_matrix.sh -s 07_bin_matrices/samples/mouse_chip_samples.tsv \
       -b data/bins/mm10_10kb_bins.bed -o data/matrices/mouse_10kb_chip_counts.txt -d CHIP/

# DNA methylation (mean per bin; % or fraction BigWigs both scaled to 0-1, see log)
sbatch 07_bin_matrices/build_methylation_matrix.sh -s 07_bin_matrices/samples/mouse_methylation_samples.tsv \
       -b data/bins/mm10_10kb_bins.bed -o data/matrices/mouse_10kb_5mC.txt -d 5mC/ -z
```

For Fig. 3 use the 5 kb bins and write `data/matrices/mouse_5kb_chip_counts.txt` and
`data/matrices/mouse_5kb_MII_5mC.txt`. Counting uses `bedtools multicov` defaults (reads overlapping a bin;
duplicate-flagged reads not counted). DNA methylation per bin is the mean of the covered CpGs
(`multiBigwigSummary BED-file`); with `-z`, bins without CpG data are 0.

## Contact
Adnan Hashim: (adnanha@uio.no)
Correspondence: Mads Lerdrup (mlerdrup@sund.ku.dk) or John Arne Dahl (j.a.dahl@medisin.uio.no)
