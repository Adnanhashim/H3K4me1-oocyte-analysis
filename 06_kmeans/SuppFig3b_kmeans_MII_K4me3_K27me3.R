# =============================================================================
# SuppFig3b_kmeans_MII_K4me3_K27me3.R
#
# Supp. Fig. 3b  Clustered heatmap (k = 5) of 10 kb bins based on MII oocyte
#                H3K4me3 and H3K27me3 ChIP-seq and P10 GO H3K4me1 picoChIP-seq
#                (two replicates each).
#
# Method
#   1. Read counts in mm10 10 kb bins without chrM/chrY (263,389 bins;
#      bedtools multicov), each dataset normalized to its top 5% signal.
#   2. Bins without signal in all six datasets are removed (and, if an input
#      count file is given, bins with no input reads).
#   3. k-means, k = 5 (nstart = 25, seed 123) on the six normalized datasets.
#   4. Clusters ordered by size (smallest first): strong H3K27me3, weaker
#      H3K27me3, strong H3K4me3, weaker H3K4me3, low signal.
#   5. Heatmap of normalized signal (0-1, white-black), bins ordered by cluster.
#
# Input (data/matrices/, shared with Fig. 2f)
#   mouse_10kb_chip_counts.txt     read counts per 10 kb bin
#   (optional) input read counts per 10 kb bin, see input_file below
# Output (figures/SuppFig3/SuppFig3b/)
#   SuppFig3b_kmeans_heatmap.tiff     the heatmap
#   kmeans_clusters.txt               bin -> cluster (heatmap order)
#   kmeans_cluster_summary.csv        cluster sizes and mean signal per dataset
#   normalized_matrix_top5.txt        top-5% normalized signal, all datasets
#   sessionInfo.txt
#
# Usage (RStudio): open H3K4me1-oocyte-analysis.Rproj and click "Source".
# Requires: pheatmap, RColorBrewer
# =============================================================================

source("06_kmeans/kmeans_clustering.R")

# ---------------------------- input / output ---------------------------------
matrix_file <- "data/matrices/mouse_10kb_chip_counts.txt"
# Optional: input (no-antibody) read counts per 10 kb bin. Bins with < 1 read in
# any input column are removed (low/no mappability), as described in the
# Methods. NULL = no input filter.
input_file  <- NULL
out_dir     <- "figures/SuppFig3/SuppFig3b"

# ---------------------------- parameters -------------------------------------
cluster_cols <- c("p10A_H3K4me1_rep1", "p10B_H3K4me1_rep2",
                  "MII_H3K4me3_rep1_GSM1871960", "MII_H3K4me3_rep2_GSM1871961",
                  "MII_H3K27me3_rep1_GSM2082667", "MII_H3K27me3_rep2_GSM2082668")
k    <- 5
seed <- 123

# ---------------------------- run --------------------------------------------
if (!file.exists(matrix_file)) {
    stop("File not found: ", matrix_file, "\nWorking directory: ", getwd(), call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(...) file.path(out_dir, paste0(...))

counts <- read_bin_matrix(matrix_file)
norm   <- normalize_by_top_fraction(counts, fraction = 0.05)
write.table(norm, out("normalized_matrix_top5.txt"), sep = "\t", quote = FALSE)

mat <- filter_bins_with_signal(norm[, cluster_cols], cluster_cols)
if (!is.null(input_file)) {
    input <- read_bin_matrix(input_file)
    has_input <- rownames(input)[apply(input >= 1, 1, all)]
    mat <- mat[rownames(mat) %in% has_input, , drop = FALSE]
}
message("Bins: ", nrow(norm), " total, ", nrow(mat), " with signal")

clusters <- run_kmeans(mat, k = k, nstart = 25, seed = seed, order_by = "size")
write.table(clusters, out("kmeans_clusters.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_tab <- data.frame(cluster = names(attr(clusters, "sizes")),
                          n_bins  = as.integer(attr(clusters, "sizes")),
                          round(attr(clusters, "centers"), 3), check.names = FALSE)
write.csv(summary_tab, out("kmeans_cluster_summary.csv"), row.names = FALSE)
print(summary_tab, row.names = FALSE)

plot_cluster_heatmap(norm[, cluster_cols], clusters,
                     file = out("SuppFig3b_kmeans_heatmap.tiff"),
                     title = "10 kb bins, top 5% normalized, k = 5")

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
