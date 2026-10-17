# =============================================================================
# Fig4b_bovine_grouped_heatmap.R
#
# Fig. 4b  Grouped heatmap of bovine oocyte H3K4me1 (GO1, GO2, FGO), FGO
#          H3K4me3 (rep1, rep2), FGO H3K27me3 (rep1, rep2) and FGO DNA
#          methylation in bosTau8 10 kb bins, grouped by FGO H3K27me3 and
#          FGO H3K4me3 status.
#
# Method
#   1. Read counts in bosTau8 10 kb bins, each dataset normalized to its top
#      5% signal (00_shared/normalization.R).
#   2. Positive/negative thresholds by mixture modelling
#      (05_mixture_model/mixture_model.R) on the mean of the two replicates:
#      three-component normal mixture on the normalized signal (the third
#      component takes the few bins with very high signal); threshold = mean
#      of the two lowest component means.
#   3. Bins grouped as: H3K27me3+ H3K4me3-  |  H3K27me3+ H3K4me3+ (overlap)  |
#      H3K27me3- H3K4me3+  |  negative. Random (seeded) order within groups.
#   4. Heatmap (0-1, white-black); FGO DNA methylation shown as fraction; the
#      binary mixture-model calls (H3K4me3, H3K27me3) as the last two columns.
#
# Input (data/Fig4/)
#   bovine_10kb_chip_counts.txt   10 kb read counts (H3K4me1, H3K4me3, H3K27me3,
#                                 H3K36me3)
#   bovine_10kb_FGO_5mC.txt      10 kb FGO DNA methylation (%)
# Output (figures/Fig4/Fig4b/)
#   Fig4b_grouped_heatmap.tiff        Fig. 4b
#   Fig4b_bins_groups_values.txt      bin, group, plotted values and binary calls
#                                     (heatmap order)
#   Fig4b_group_sizes.csv             number of bins per group
#   thresholds.csv                    mixture-model thresholds
#   mixture_model_parameters.csv      mixture-model parameters
#   mixture_model_fits.pdf            fitted mixture per dataset
#   sessionInfo.txt
#
# Usage (RStudio): open H3K4me1-oocyte-analysis.Rproj and click "Source".
# Requires: pheatmap, RColorBrewer, mixtools, fBasics
# =============================================================================

suppressPackageStartupMessages({
    library(pheatmap)
    library(RColorBrewer)
})
for (f in c("00_shared/normalization.R", "05_mixture_model/mixture_model.R")) {
    if (!file.exists(f)) {
        stop("Shared script not found: ", f, "\nWorking directory: ", getwd(),
             "\nOpen H3K4me1-oocyte-analysis.Rproj (repository root) before running.", call. = FALSE)
    }
    source(f)
}
if (!"low_mid_mean" %in% eval(formals(fit_mixture_models)$threshold_rule)) {
    stop("05_mixture_model/mixture_model.R is an older version - replace it with the current one.",
         call. = FALSE)
}

# ---------------------------- input / output ---------------------------------
data_dir    <- "data/Fig4"
counts_file <- file.path(data_dir, "bovine_10kb_chip_counts.txt")
meth_file   <- file.path(data_dir, "bovine_10kb_FGO_5mC.txt")
out_dir     <- "figures/Fig4/Fig4b"
seed        <- 123

# Column names in the count matrix -> labels used in the figure
chip_cols <- c(GO1_H3K4me1        = "NSN_K4me1_rep1",
               GO2_H3K4me1        = "NSN_K4me1_rep2",
               FGO_H3K4me1        = "SN_K4me1_rep1",
               FGO_H3K4me3_rep1   = "FGO_H3K4me3_rep1",
               FGO_H3K4me3_rep2   = "FGO_H3K4me3_rep2",
               FGO_H3K27me3_rep1  = "FGO_H3K27me3_rep1",
               FGO_H3K27me3_rep2  = "FGO_H3K27me3_rep2")

# ---------------------------- checks -----------------------------------------
for (f in c(counts_file, meth_file)) {
    if (!file.exists(f)) stop("File not found: ", f, "\nWorking directory: ", getwd(), call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(...) file.path(out_dir, paste0(...))

# ---------------------------- 1. load and normalize --------------------------
counts <- read_bin_matrix(counts_file)
missing_cols <- setdiff(chip_cols, colnames(counts))
if (length(missing_cols) > 0) {
    stop("Columns not found in ", counts_file, ": ", paste(missing_cols, collapse = ", "), call. = FALSE)
}
norm <- normalize_by_top_fraction(counts, fraction = 0.05)   # all datasets, as in the original analysis
bins <- rownames(norm)

meth <- read_bin_matrix(meth_file)
values <- norm[, chip_cols]
colnames(values) <- names(chip_cols)
values$FGO_DNAme <- meth[match(bins, rownames(meth)), 1] / 100   # % -> fraction

# ---------------------------- 2. mixture modelling ---------------------------
mix_input <- data.frame(
    FGO_H3K4me3_avg  = (values$FGO_H3K4me3_rep1 + values$FGO_H3K4me3_rep2) / 2,
    FGO_H3K27me3_avg = (values$FGO_H3K27me3_rep1 + values$FGO_H3K27me3_rep2) / 2,
    row.names = bins
)
mix <- fit_mixture_models(mix_input, k = 3, seed = seed,
                          threshold_rule = "low_mid_mean", scale = "linear",
                          start = "quantile", start_sd = c(0.05, 0.1, 0.2, 0.5),
                          plot_file = out("mixture_model_fits.pdf"),
                          density_adjust = 0.1, plot_xlim = c(0, 5))
write.csv(mix$table, out("mixture_model_parameters.csv"), row.names = FALSE)
thr <- setNames(exp(mix$table$threshold_ln), mix$table$dataset)    # normalized scale
write.csv(data.frame(dataset = names(thr), threshold = thr), out("thresholds.csv"), row.names = FALSE)
print(round(thr, 4))

# ---------------------------- 3. groups --------------------------------------
k4me3  <- mix_input$FGO_H3K4me3_avg  > thr["FGO_H3K4me3_avg"]
k27me3 <- mix_input$FGO_H3K27me3_avg > thr["FGO_H3K27me3_avg"]
group  <- ifelse(k27me3 & !k4me3, "H3K27me3",
          ifelse(k27me3 &  k4me3, "Overlap",
          ifelse(!k27me3 & k4me3, "H3K4me3", "Negative")))
group_levels <- c("H3K27me3", "Overlap", "H3K4me3", "Negative")

set.seed(seed)
shuffled <- sample(seq_along(bins))
o <- shuffled[order(match(group[shuffled], group_levels))]

sizes <- as.data.frame(table(factor(group, levels = group_levels)))
colnames(sizes) <- c("group", "n_bins")
write.csv(sizes, out("Fig4b_group_sizes.csv"), row.names = FALSE)
print(sizes, row.names = FALSE)

# ---------------------------- 4. heatmap -------------------------------------
group_labels <- setNames(paste0(sizes$group, " (n = ", formatC(sizes$n_bins, format = "d", big.mark = ","), ")"),
                         as.character(sizes$group))
values$Binary_FGO_H3K4me3  <- as.integer(k4me3)
values$Binary_FGO_H3K27me3 <- as.integer(k27me3)
hm <- as.matrix(values[o, ])
write.table(data.frame(bin = bins[o], group = group[o], hm, check.names = FALSE),
            out("Fig4b_bins_groups_values.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

pheatmap(hm, cluster_cols = FALSE, cluster_rows = FALSE,
         color = colorRampPalette(brewer.pal(9, "Greys"))(100),
         breaks = seq(0, 1, length.out = 101),
         annotation_row = data.frame(group = factor(group_labels[group[o]], levels = group_labels),
                                     row.names = rownames(hm)),
         show_rownames = FALSE, gaps_col = c(3, 5, 7, 8),
         main = "bosTau8 10 kb bins, top 5% normalized",
         filename = out("Fig4b_grouped_heatmap.tiff"), width = 6, height = 10)

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
