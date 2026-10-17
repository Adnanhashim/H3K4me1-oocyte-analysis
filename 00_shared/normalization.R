# =============================================================================
# normalization.R
#
# Normalization of signal in genomic bins (bins x datasets), shared by the
# mixture-modelling (05_mixture_model) and k-means (06_kmeans) scripts.
# =============================================================================

#' Read a bin x dataset matrix (first column = bin name "chr_start_end").
#' Uses data.table::fread when installed (much faster for 400 bp matrices).
read_bin_matrix <- function(file) {
    if (requireNamespace("data.table", quietly = TRUE)) {
        d <- as.data.frame(data.table::fread(file, sep = "\t", header = TRUE, check.names = FALSE))
        m <- d[, -1, drop = FALSE]
        rownames(m) <- d[[1]]
    } else {
        m <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE,
                        row.names = 1, comment.char = "")
    }
    colnames(m) <- trimws(colnames(m))
    m[] <- lapply(m, as.numeric)
    m
}

#' Divide each column by its 95th percentile (Fig. 1e/1f).
normalize_by_quantile <- function(mat, probs = 0.95) {
    q <- apply(mat, 2, function(x) quantile(x, probs = probs))
    sweep(mat, 2, q, "/")
}

#' Divide each column by the value at rank fraction * n when sorted
#' decreasingly ("top 5% normalization"; k-means and grouped heatmaps).
normalize_by_top_fraction <- function(mat, fraction = 0.05) {
    index <- fraction * nrow(mat)
    out <- as.data.frame(apply(mat, 2, function(col) {
        col / sort(col, decreasing = TRUE)[index]
    }))
    rownames(out) <- rownames(mat)
    out
}

#' Natural log per column; bins without signal (log = -Inf) get zero_value.
log_transform <- function(mat, zero_value = -Inf) {
    out <- as.data.frame(lapply(mat, function(x) {
        lx <- log(x)
        lx[!is.finite(lx)] <- zero_value
        lx
    }))
    rownames(out) <- rownames(mat)
    out
}
