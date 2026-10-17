# =============================================================================
# Fig4e_human_subclustering_400bp.R
#
# Fig. 4e and Supp. Fig. 4c
# Fig. 4e  Subclustering of the H3K4me3-positive / H3K27me3-negative 10 kb bins
#          from Fig. 4d at 400 bp resolution: MI/II H3K4me1, 4-cell H3K4me3
#          (rep1, rep2) and MII DNA methylation.
#
# Method
#   1. Read counts in hg38 400 bp bins, each dataset normalized to its top 5%
#      signal genome-wide (00_shared/normalization.R).
#   2. Each H3K4me3+/H3K27me3- 10 kb bin (Fig. 4d) is split into 25 x 400 bp
#      bins; for every dataset the 25 values form one row segment
#      (4 datasets x 25 = 100 values per 10 kb bin).
#   3. k-means, k = 10 (nstart = 25, seed 123; 06_kmeans/kmeans_clustering.R).
#   4. Clusters labelled by their mean profile and ordered as in the
#      manuscript: broad H3K4me3 (intermediate, high), focal H3K4me3 peaks
#      (3' to 5'), high DNA methylation, low signal, rep2-only, outliers.
#   5. Heatmap (0-1, white-black), 10 kb bins ordered by cluster.
#
# Input
#   figures/Fig4/Fig4d/K4me3pos_K27me3neg_bins.txt       from Fig4d script
#   data/Fig4/human_400bp_chip_counts.txt   400 bp read counts
#   data/Fig4/human_400bp_MII_5mC.txt              400 bp DNA methylation (fraction)
# Output (figures/Fig4/Fig4e/)
#   Fig4e_subclustering_heatmap.tiff  Fig. 4e
#   kmeans_clusters.txt               10 kb bin -> cluster (heatmap order)
#   kmeans_cluster_summary.csv        cluster sizes and mean signal per dataset
#   sessionInfo.txt
# Output (figures/SuppFig4/SuppFig4c/)
#   SuppFig4c_subclustering_heatmap.tiff  same clusters with GO1/GO2 H3K4me1 and
#                                         GV H3K27me3 added
#
# Usage (RStudio): open H3K4me1-oocyte-analysis.Rproj and click "Source".
# Requires: pheatmap, RColorBrewer
# =============================================================================

source("06_kmeans/kmeans_clustering.R")

# ---------------------------- input / output ---------------------------------
bins_file   <- "figures/Fig4/Fig4d/K4me3pos_K27me3neg_bins.txt"
counts_file <- "data/Fig4/human_400bp_chip_counts.txt"
meth_file   <- "data/Fig4/human_400bp_MII_5mC.txt"
out_dir     <- "figures/Fig4/Fig4e"

# ---------------------------- parameters -------------------------------------
# Column names in the 400 bp count matrix -> labels used in the figure
chip_cols <- c(MI_II_H3K4me1   = "10-hMA-K4m1-22_S8",
               C4_H3K4me3_rep1 = "hs_4cell_2PN_H3K4me3_rep1",
               C4_H3K4me3_rep2 = "hs_4cell_2PN_H3K4me3_rep2")
# Additional datasets shown only in Supp. Fig. 4c (same clusters as Fig. 4e)
supp_cols <- c(hPriA = "13-hPriA-K4m1", hPriB = "14-hPriB-K4m1",
               hSecA = "15-hSecA-K4m1", hSecB = "16-hSecB-K4m1",
               GV_H3K27me3 = "hs_GV_H3K27me3_rep2")
supp_out_dir <- "figures/SuppFig4/SuppFig4c"
sub_bin <- 400
bin     <- 10000
k       <- 10
seed    <- 123

# ---------------------------- checks -----------------------------------------
for (f in c(bins_file, counts_file, meth_file)) {
    if (!file.exists(f)) stop("File not found: ", f, "\nWorking directory: ", getwd(), call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(...) file.path(out_dir, paste0(...))

# ---------------------------- 1. 400 bp signal -------------------------------
counts <- read_bin_matrix(counts_file)
missing_cols <- setdiff(c(chip_cols, supp_cols), colnames(counts))
if (length(missing_cols) > 0) {
    stop("Columns not found in ", counts_file, ": ", paste(missing_cols, collapse = ", "),
         "\nColumns present: ", paste(colnames(counts), collapse = ", "), call. = FALSE)
}
chip <- normalize_by_top_fraction(counts[, chip_cols, drop = FALSE], fraction = 0.05)
colnames(chip) <- names(chip_cols)

meth <- read_bin_matrix(meth_file)
signal_400 <- cbind(chip, MII_DNAme = meth[match(rownames(chip), rownames(meth)), 1])
datasets <- colnames(signal_400)

# ---------------------------- 2. 10 kb bins x 25 positions -------------------
bins10 <- readLines(bins_file)
parts  <- do.call(rbind, strsplit(bins10, "_"))
start  <- as.numeric(parts[, 2]); end <- as.numeric(parts[, 3])
bins10 <- bins10[end - start == bin]                    # full-length bins only
parts  <- parts[end - start == bin, , drop = FALSE]

n_sub  <- bin / sub_bin
offset <- (seq_len(n_sub) - 1) * sub_bin
sub_names <- outer(seq_along(bins10), seq_len(n_sub), function(i, j) {
    s <- as.numeric(parts[i, 2]) + offset[j]
    paste(parts[i, 1], s, s + sub_bin, sep = "_")
})
row_idx <- matrix(match(sub_names, rownames(signal_400)), nrow = length(bins10))

mat <- do.call(cbind, lapply(datasets, function(d) {
    m <- matrix(signal_400[[d]][row_idx], nrow = length(bins10))
    colnames(m) <- paste0(d, "_", offset); m
}))
rownames(mat) <- bins10
complete <- complete.cases(mat)
mat <- mat[complete, , drop = FALSE]
message("10 kb bins: ", length(bins10), " used: ", nrow(mat),
        " (", sum(!complete), " with missing 400 bp values removed)")

# ---------------------------- 3-4. k-means -----------------------------------
# Clusters are labelled by their mean profile and placed in the manuscript
# order (k-means numbers are arbitrary):
#   broad_medium  broad 4-cell H3K4me3, intermediate level
#   broad_high    broad 4-cell H3K4me3, high level
#   peak          focal 4-cell H3K4me3 peak; ordered by peak position,
#                 right (3') to left (5') within the 10 kb bin
#   high_DNAme    low H3K4me3, high MII DNA methylation
#   low           low H3K4me3 and DNA methylation
#   rep2_only     signal only in 4-cell H3K4me3 rep2
#   outlier       cluster with < 10 bins
label_order <- c("broad_medium", "broad_high", "peak", "high_DNAme", "low", "rep2_only", "outlier")

cluster_profile <- function(centers, sizes) {
    prof <- function(d) centers[, grep(paste0("^", d, "_"), colnames(centers)), drop = FALSE]
    rep1 <- prof("C4_H3K4me3_rep1"); rep2 <- prof("C4_H3K4me3_rep2")
    k4me3 <- (rep1 + rep2) / 2
    data.frame(
        k4me3_mean = rowMeans(k4me3),
        sharpness  = apply(k4me3, 1, max) / rowMeans(k4me3),   # max / mean of profile
        peak_pos   = apply(k4me3, 1, which.max),               # 1..25, 5' -> 3'
        rep1_rep2  = rowMeans(rep1) / pmax(rowMeans(rep2), 1e-9),
        dname_mean = rowMeans(prof("MII_DNAme")),
        size       = sizes
    )
}

label_clusters <- function(p) {
    ifelse(p$size < 10, "outlier",
    ifelse(p$rep1_rep2 < 0.2, "rep2_only",
    ifelse(p$sharpness >= 1.2, "peak",
    ifelse(p$k4me3_mean >= 0.9, "broad_high",
    ifelse(p$dname_mean >= 0.5, "high_DNAme",
    ifelse(p$k4me3_mean >= 0.55, "broad_medium", "low"))))))
}

order_like_manuscript <- function(centers, sizes) {
    p <- cluster_profile(centers, sizes)
    lab <- label_clusters(p)
    order(match(lab, label_order), -p$peak_pos, -p$k4me3_mean)
}

clusters <- run_kmeans(as.data.frame(mat), k = k, nstart = 25, seed = seed,
                       order_by = order_like_manuscript)
write.table(clusters, out("kmeans_clusters.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

centers <- attr(clusters, "centers")
profile <- cluster_profile(centers, attr(clusters, "sizes"))
summary_tab <- data.frame(cluster = rownames(centers),
                          pattern = label_clusters(profile),
                          k4me3_peak_position = profile$peak_pos,
                          n_bins = as.integer(attr(clusters, "sizes")),
                          sapply(datasets, function(d) round(rowMeans(centers[, grep(paste0("^", d, "_"), colnames(centers))]), 3)),
                          check.names = FALSE)
write.csv(summary_tab, out("kmeans_cluster_summary.csv"), row.names = FALSE)
print(summary_tab, row.names = FALSE)

# ---------------------------- 5. heatmap -------------------------------------
# One label per dataset, under the middle of its 25 x 400 bp block
dataset_labels <- c(MI_II_H3K4me1 = "H3K4me1\nMI/II", C4_H3K4me3_rep1 = "H3K4me3\n4C rep1",
                    C4_H3K4me3_rep2 = "H3K4me3\n4C rep2", MII_DNAme = "DNAme\nMII")
block_labels <- rep("", ncol(mat))
block_labels[(seq_along(datasets) - 1) * n_sub + ceiling(n_sub / 2)] <- dataset_labels[datasets]

plot_cluster_heatmap(as.data.frame(mat), clusters,
                     file = out("Fig4e_subclustering_heatmap.tiff"),
                     title = "H3K4me3+/H3K27me3- 10 kb bins, 400 bp resolution",
                     gaps_col = n_sub * seq_len(length(datasets) - 1),
                     labels_col = block_labels, angle_col = "0",
                     cluster_names = paste0(rownames(centers), " ", summary_tab$pattern))

# ---------------------------- Supp. Fig. 4c ----------------------------------
# Same 10 kb bins, clusters and order as Fig. 4e, with GO1/GO2 H3K4me1 (mean of
# two oocytes each) and GV H3K27me3 added. Clustering is not repeated.
supp <- normalize_by_top_fraction(counts[, supp_cols, drop = FALSE], fraction = 0.05)
colnames(supp) <- names(supp_cols)
signal_supp <- data.frame(
    GO1_H3K4me1     = (supp$hPriA + supp$hPriB) / 2,
    GO2_H3K4me1     = (supp$hSecA + supp$hSecB) / 2,
    MI_II_H3K4me1   = signal_400$MI_II_H3K4me1,
    C4_H3K4me3_rep1 = signal_400$C4_H3K4me3_rep1,
    C4_H3K4me3_rep2 = signal_400$C4_H3K4me3_rep2,
    GV_H3K27me3     = supp$GV_H3K27me3,
    MII_DNAme       = signal_400$MII_DNAme,
    row.names = rownames(signal_400)
)
supp_datasets <- colnames(signal_supp)
idx_used <- row_idx[complete, , drop = FALSE]
mat_supp <- do.call(cbind, lapply(supp_datasets, function(d) {
    m <- matrix(signal_supp[[d]][idx_used], nrow = nrow(idx_used))
    colnames(m) <- paste0(d, "_", offset); m
}))
rownames(mat_supp) <- rownames(mat)

supp_labels <- c(GO1_H3K4me1 = "H3K4me1\nGO1", GO2_H3K4me1 = "H3K4me1\nGO2",
                 dataset_labels, GV_H3K27me3 = "H3K27me3\nGV")
supp_block_labels <- rep("", ncol(mat_supp))
supp_block_labels[(seq_along(supp_datasets) - 1) * n_sub + ceiling(n_sub / 2)] <- supp_labels[supp_datasets]

dir.create(supp_out_dir, recursive = TRUE, showWarnings = FALSE)
plot_cluster_heatmap(as.data.frame(mat_supp), clusters,
                     file = file.path(supp_out_dir, "SuppFig4c_subclustering_heatmap.tiff"),
                     title = "H3K4me3+/H3K27me3- 10 kb bins, 400 bp resolution",
                     gaps_col = n_sub * seq_len(length(supp_datasets) - 1),
                     labels_col = supp_block_labels, angle_col = "0",
                     cluster_names = paste0(rownames(centers), " ", summary_tab$pattern),
                     width = 9)

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
