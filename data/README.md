# data/

Input matrices for the figure scripts (not tracked in git). Raw data: GEO (see the main README and
Supp. Table 1). The mouse bin matrices are made on the server with `07_bin_matrices/`; one matrix per bin
size holds all samples and is shared between figures.

```
data/
├── matrices/  mouse_10kb_chip_counts.txt   mm10 10 kb read counts, all ChIP + RNA samples   (Fig. 2f, Supp. Fig. 3b)
│              mouse_10kb_5mC.txt           mm10 10 kb DNA methylation, fraction             (Fig. 2f)
│              mouse_5kb_chip_counts.txt    mm10 5 kb read counts, all ChIP samples          (Fig. 3a-c)
│              mouse_5kb_MII_5mC.txt        mm10 5 kb MII DNA methylation, fraction          (Fig. 3b)
├── Fig1/      Fig1F_mm10_10KB_bins_scatterplot.txt            10 kb read counts, H3K4me1/H3K27ac (Fig. 1e/1f)
│              Fig1F_mm10_10KB_bins_TSS_distalpeaks_input.txt  TSS-distal 10 kb bins
└── Fig4/      human_10kb_*.txt, human_400bp_*.txt             hg38 matrices (Fig. 4d, 4e, Supp. Fig. 4c)
               bovine_10kb_*.txt                               bosTau8 matrices (Fig. 4b)
```
