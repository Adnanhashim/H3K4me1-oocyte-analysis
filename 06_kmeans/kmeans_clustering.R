# =============================================================================
# kmeans_clustering.R
#
# Shared functions for k-means clustering of genomic bins and the clustered
# heatmaps. All k-means clustering in the manuscript uses these functions;
# each figure script sources this file and sets its own parameters:
#
#   Figure     Bins    Normalization   Clustered on                       k
#   Supp. 3b   10 kb   top 5%          P10 GO H3K4me1, MII H3K4me3,       5
#                                      MII H3K27me3 (2 replicates each)
#   4e         400 bp  top 5%          human MI/II H3K4me1, 4-cell        10
#              within                  H3K4me3 (rep1, rep2), MII DNA
#              10 kb                   methylation; 25 x 400 bp per
#                                      H3K4me3+/H3K27me3- 10 kb bin (Fig. 4d)
#   (further figures are added here as they are converted)
#
# Method:
#   1. read counts per bin (bedtools multicov) normalized per dataset to the
#      top 5% signal (00_shared/normalization.R)
#   2. bins without signal in any clustering dataset are removed
#   3. k-means (stats::kmeans, nstart = 25, fixed seed)
#   4. clusters ordered by a fixed rule (default: cluster size) so the order
#      does not depend on the arbitrary cluster numbers
#   5. heatmap (pheatmap, white-black, 0-1), bins ordered by cluster
#
# Usage: source("06_kmeans/kmeans_clustering.R") from a figure script.
# Requires: pheatmap, RColorBrewer
# =============================================================================

suppressPackageStartupMessages({
    library(pheatmap)
    library(RColorBrewer)
})
source("00_shared/normalization.R")

#' Keep bins with signal (> 0) in at least one of the given columns.
filter_bins_with_signal <- function(mat, cols = colnames(mat)) {
    mat[apply(mat[, cols, drop = FALSE] > 0, 1, any), , drop = FALSE]
}

#' k-means with a fixed seed; clusters renumbered 1..k by 'order_by'.
#'
#' @param mat       bins x datasets (clustering columns only)
#' @param order_by  "size" (smallest cluster first), a function(centers, sizes)
#'                  returning the cluster order, or an integer vector giving the
#'                  original cluster numbers in the wanted order
#' @return data.frame(bin, cluster) in heatmap order, plus attributes
#'         'centers' and 'sizes'
run_kmeans <- function(mat, k, nstart = 25, seed = 123, order_by = "size") {
    set.seed(seed)
    km <- kmeans(mat, centers = k, nstart = nstart, iter.max = 100)

    new_order <- if (identical(order_by, "size")) {
        order(km$size)
    } else if (is.function(order_by)) {
        order_by(km$centers, km$size)
    } else order_by
    relabel <- setNames(seq_along(new_order), new_order)
    cl <- relabel[as.character(km$cluster)]

    # bins ordered by cluster; random (seeded) order within a cluster
    set.seed(seed)
    shuffled <- sample(seq_len(nrow(mat)))
    o <- shuffled[order(cl[shuffled])]

    res <- data.frame(bin = rownames(mat)[o], cluster = unname(cl[o]),
                      stringsAsFactors = FALSE)
    centers <- km$centers[new_order, , drop = FALSE]
    rownames(centers) <- paste0("Cluster", seq_len(k))
    attr(res, "centers") <- centers
    attr(res, "sizes") <- setNames(km$size[new_order], rownames(centers))
    res
}

#' Clustered heatmap (pheatmap) of selected columns for bins in cluster order.
#'
#' @param values    bins x datasets to display (row names = bins)
#' @param clusters  output of run_kmeans()
#' @param gaps_col  column positions for vertical gaps
#' @param row_gaps  TRUE = white gap between clusters
#' @param labels_col  optional column labels (e.g. one label per dataset block)
#' @param cluster_names  optional names per cluster (in cluster order); the
#'                  annotation shows "<name> (n = <bins>)"
plot_cluster_heatmap <- function(values, clusters, file, title = "",
                                 gaps_col = NULL, row_gaps = FALSE,
                                 show_colnames = TRUE, cluster_names = NULL,
                                 labels_col = NULL, angle_col = "270",
                                 width = 6, height = 10) {
    hm <- as.matrix(values[match(clusters$bin, rownames(values)), , drop = FALSE])
    rownames(hm) <- clusters$bin
    # Annotation label per cluster, with the number of bins: "Cluster1 (n = 6,412)"
    ids <- sort(unique(clusters$cluster))
    n <- table(clusters$cluster)[as.character(ids)]
    base <- if (is.null(cluster_names)) paste0("Cluster", ids) else cluster_names[ids]
    labels <- paste0(base, " (n = ", formatC(as.integer(n), format = "d", big.mark = ","), ")")
    ann <- data.frame(cluster = factor(labels[match(clusters$cluster, ids)], levels = labels),
                      row.names = clusters$bin)
    gaps_row <- cumsum(table(clusters$cluster))
    gaps_row <- if (row_gaps) gaps_row[-length(gaps_row)] else NULL

    pheatmap(hm, cluster_cols = FALSE, cluster_rows = FALSE,
             color = colorRampPalette(brewer.pal(9, "Greys"))(100),
             breaks = seq(0, 1, length.out = 101),
             annotation_row = ann, show_rownames = FALSE, show_colnames = show_colnames,
             labels_col = if (is.null(labels_col)) colnames(hm) else labels_col,
             angle_col = angle_col,
             gaps_row = gaps_row, gaps_col = gaps_col,
             main = title, filename = file, width = width, height = height)
}
