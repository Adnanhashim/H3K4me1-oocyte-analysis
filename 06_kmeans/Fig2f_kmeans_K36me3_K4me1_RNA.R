# =============================================================================
# Fig2f_kmeans_K36me3_K4me1_RNA.R
#
# Fig. 2f  Clustered heatmap (k = 15) of mm10 10 kb bins based on P10 GO,
#          P14 GO and FGO H3K36me3 (rep1, rep2) and RNA-seq, and P10 GO and
#          4-week NSN H3K4me1 picoChIP-seq (rep1, rep2). DNA methylation in
#          P10 GO, P15 GO and FGO is displayed on the resulting clusters.
#
# Method
#   1. Read counts in mm10 10 kb bins without chrM/chrY (263,389 bins;
#      07_bin_matrices/), each dataset normalized to its top 5% signal.
#   2. Bins without signal in all clustering datasets are removed.
#   3. k-means, k = 15 (nstart = 25, seed 123) on H3K36me3, H3K4me1 and RNA.
#   4. Clusters ordered as in the manuscript: H3K36me3-marked clusters (mean
#      H3K36me3 > k36_split) first, ordered by RNA (highest first); then the
#      H3K36me3-low clusters, ordered by H3K4me1 (highest first).
#   5. Heatmap (0-1, white-black): H3K36me3 | H3K4me1 | RNA | DNA methylation
#      (fraction; not used for clustering).
#
# Input (data/matrices/, shared with Supp. Fig. 3b), made with 07_bin_matrices/ (column names = sample sheet labels)
#   mouse_10kb_chip_counts.txt   read counts per 10 kb bin (ChIP and RNA)
#   mouse_10kb_5mC.txt           DNA methylation (fraction) per 10 kb bin
# Output (figures/Fig2/Fig2f/)
#   Fig2f_kmeans_heatmap.tiff    the heatmap
#   kmeans_clusters.txt          bin -> cluster (heatmap order)
#   kmeans_cluster_summary.csv   cluster sizes and mean signal per dataset
#                                (incl. DNA methylation)
#   sessionInfo.txt
#
# Usage (RStudio): open H3K4me1-oocyte-analysis.Rproj and click "Source".
# Requires: pheatmap, RColorBrewer
# =============================================================================

source("06_kmeans/kmeans_clustering.R")

# ---------------------------- input / output ---------------------------------
chip_file <- "data/matrices/mouse_10kb_chip_counts.txt"
meth_file <- "data/matrices/mouse_10kb_5mC.txt"
out_dir   <- "figures/Fig2/Fig2f"
k         <- 15
seed      <- 123
k36_split <- 0.3     # mean normalized H3K36me3 separating marked / low clusters

# ---------------------------- columns (matrix labels) ------------------------
# Plot name = column in the matrix (labels from 07_bin_matrices/samples/).
# GSM numbers from Supp. Table 1. Check the P10 H3K4me1 samples: the Results state that BL6 GO H3K4me1 data
# were used here (the CD1 p10A/p10B labels are set as default).
k36_cols <- c(GO_P10_rep1 = "P10_H3K36me3_rep1_GSM3084618",
              GO_P10_rep2 = "P10_H3K36me3_rep2_GSM3084619",
              GO_P14_rep1 = "P14_H3K36me3_rep1_GSM3084620",
              GO_P14_rep2 = "P14_H3K36me3_rep2_GSM3084621",
              FGO_rep1    = "FGO_H3K36me3_rep1_GSM3084622",
              FGO_rep2    = "FGO_H3K36me3_rep2_GSM3084623")
k4me1_cols <- c(GO_P10_rep1 = "p10A_H3K4me1_rep1",
                GO_P10_rep2 = "p10B_H3K4me1_rep2",
                NSN_4w_rep1 = "3-NSNf-m1-50_S23",
                NSN_4w_rep2 = "NSNB_m1-50-40")
rna_cols  <- c(GO_P10 = "P10_RNA_GSM1845292",
               GO_P14 = "P14_RNA_GSM1845293",
               FGO    = "FGO_RNA_GSM1845294")
# DNA methylation: one or more columns per stage (shown as their mean)
meth_cols <- list(GO_P10 = "P10_5mC_GSM2690945_46",
                  GO_P15 = "P15_5mC_GSE72784",
                  FGO    = "FGO_5mC_GSM4832691_92_93")

# ---------------------------- checks -----------------------------------------
for (f in c(chip_file, meth_file)) {
    if (!file.exists(f)) stop("File not found: ", f, "\nWorking directory: ", getwd(), call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(...) file.path(out_dir, paste0(...))

counts  <- read_bin_matrix(chip_file)
meth    <- read_bin_matrix(meth_file)
cluster_cols <- c(k36_cols, k4me1_cols, rna_cols)
missing <- c(setdiff(cluster_cols, colnames(counts)), setdiff(unlist(meth_cols), colnames(meth)))
if (length(missing) > 0) {
    stop("Columns not found in the input matrices: ", paste(missing, collapse = ", "),
         "\nEdit the column names at the top of this script.", call. = FALSE)
}

# ---------------------------- 1. normalize -----------------------------------
norm <- normalize_by_top_fraction(counts[, cluster_cols], fraction = 0.05)
colnames(norm) <- c(paste0("H3K36me3_", names(k36_cols)),
                    paste0("H3K4me1_",  names(k4me1_cols)),
                    paste0("RNA_",      names(rna_cols)))
meth <- meth[match(rownames(norm), rownames(meth)), unlist(meth_cols), drop = FALSE]
meth[is.na(meth)] <- 0
meth[] <- lapply(meth, function(x) if (max(x) > 1) x / 100 else x)   # % -> fraction, per file
meth <- as.data.frame(lapply(meth_cols, function(cols) rowMeans(meth[, cols, drop = FALSE])),
                      row.names = rownames(norm))
colnames(meth) <- paste0("DNAme_", names(meth_cols))

# ---------------------------- 2. filter + k-means ----------------------------
mat <- filter_bins_with_signal(norm)
message("Bins: ", nrow(norm), " total, ", nrow(mat), " with signal")

k36_names <- paste0("H3K36me3_", names(k36_cols))
k4_names  <- paste0("H3K4me1_",  names(k4me1_cols))
rna_names <- paste0("RNA_",      names(rna_cols))
manuscript_order <- function(centers, sizes) {
    k36 <- rowMeans(centers[, k36_names, drop = FALSE])
    k4  <- rowMeans(centers[, k4_names,  drop = FALSE])
    rna <- rowMeans(centers[, rna_names, drop = FALSE])
    marked <- k36 > k36_split
    c(which(marked)[order(-rna[marked])], which(!marked)[order(-k4[!marked])])
}
clusters <- run_kmeans(mat, k = k, nstart = 25, seed = seed, order_by = manuscript_order)
write.table(clusters, out("kmeans_clusters.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

values <- cbind(norm, meth)
mean_meth <- sapply(split(values[clusters$bin, colnames(meth), drop = FALSE], clusters$cluster),
                    colMeans)
summary_tab <- data.frame(cluster = names(attr(clusters, "sizes")),
                          n_bins  = as.integer(attr(clusters, "sizes")),
                          round(attr(clusters, "centers"), 3),
                          round(t(mean_meth), 3),
                          check.names = FALSE)
write.csv(summary_tab, out("kmeans_cluster_summary.csv"), row.names = FALSE)
print(summary_tab[, 1:2], row.names = FALSE)

# ---------------------------- 3. heatmap -------------------------------------
n36 <- length(k36_cols); n1 <- length(k4me1_cols); nr <- length(rna_cols)
plot_cluster_heatmap(values, clusters,
                     file = out("Fig2f_kmeans_heatmap.tiff"),
                     title = "mm10 10 kb bins, top 5% normalized, k = 15",
                     gaps_col = cumsum(c(n36, n1, nr)),
                     labels_col = colnames(values), width = 9, height = 10)

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
