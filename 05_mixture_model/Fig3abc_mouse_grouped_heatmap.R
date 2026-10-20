# =============================================================================
# Fig3abc_mouse_grouped_heatmap.R
#
# Fig. 3a  Grouped heatmap of mm10 5 kb bins based on MII oocyte H3K4me3 and
#          H3K27me3 ChIP-seq (rep1, rep2).
# Fig. 3b  MII H3K36me3 and MII DNA methylation at the groups defined in a.
# Fig. 3c  H3K4me1 (picoChIP-seq) in P10 GO, 4-week NSN, 4-week SN and MII
#          oocytes (rep1, rep2) at the groups defined in a.
#
# Method
#   1. Read counts in mm10 5 kb bins without chrM/chrY (526,765 bins;
#      07_bin_matrices/), each dataset normalized to its top 5% signal
#      (00_shared/normalization.R). For grouping, MII H3K4me3 and MII H3K27me3
#      = mean of the two replicates; MII H3K36me3 is shown as the mean of the
#      two replicates.
#   2. Positive/negative thresholds by mixture modelling
#      (05_mixture_model/mixture_model.R): two-component normal mixture on the
#      ln of the normalized signal (bins with signal); threshold = mean of the
#      two component means. The figure uses the fixed thresholds set below
#      ('fixed_thresholds'); thresholds.csv reports both.
#   3. Bins grouped as: H3K27me3+ H3K4me3-  |  H3K27me3+ H3K4me3+ (overlap)  |
#      H3K27me3- H3K4me3+  |  negative. Random (seeded) order within groups.
#   4. Three heatmaps with the same bin order (0-1, white-black); DNA
#      methylation shown as fraction.
#
# Input (data/matrices/), made with 07_bin_matrices/ (column names = sample sheet labels)
#   mouse_5kb_chip_counts.txt      read counts per 5 kb bin
#   mouse_5kb_MII_5mC.txt          DNA methylation (fraction) per 5 kb bin
# Output (figures/Fig3/)
#   Fig3a_grouped_heatmap.tiff, Fig3b_..., Fig3c_...   heatmaps
#   Fig3_bins_groups_values.txt    bin, group, plotted values (heatmap order)
#   Fig3_group_sizes.csv           number of bins per group
#   thresholds.csv                 thresholds used and mixture-model fit
#   mixture_model_parameters.csv   mixture-model parameters
#   mixture_model_fits.pdf         fitted mixture per dataset
#   sessionInfo.txt
#
# The same script runs on 10 kb bins: set bin_size <- "10kb".
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
if (!"start_sd" %in% names(formals(fit_mixture_models))) {
    stop("05_mixture_model/mixture_model.R is an older version - replace it with the current one.",
         call. = FALSE)
}

# ---------------------------- input / output ---------------------------------
bin_size  <- "5kb"
chip_file <- file.path("data/matrices", paste0("mouse_", bin_size, "_chip_counts.txt"))
meth_file <- file.path("data/matrices", paste0("mouse_", bin_size, "_MII_5mC.txt"))
out_dir   <- "figures/Fig3"
seed      <- 123


fixed_thresholds <- c(MII_H3K4me3 = 0.488, MII_H3K27me3 = 0.042)

# ---------------------------- columns (matrix labels) ------------------------
# Plot name = column in the count matrix
fig3a_cols <- c(MII_H3K4me3_rep1  = "MII_H3K4me3_rep1_GSM1871960",
                MII_H3K4me3_rep2  = "MII_H3K4me3_rep2_GSM1871961",
                MII_H3K27me3_rep1 = "MII_H3K27me3_rep1_GSM2082667",
                MII_H3K27me3_rep2 = "MII_H3K27me3_rep2_GSM2082668")
k36_cols   <- c("MII_H3K36me3_rep1_GSM3084624", "MII_H3K36me3_rep2_GSM3084625")
meth_col   <- "MII_5mC"
fig3c_cols <- c(GO_P10_rep1 = "p10A_H3K4me1_rep1",
                GO_P10_rep2 = "p10B_H3K4me1_rep2",
                NSN_4w_rep1 = "3-NSNf-m1-50_S23",
                NSN_4w_rep2 = "NSNB_m1-50-40",
                SN_4w_rep1  = "3-SN-m1-90_S3",
                SN_4w_rep2  = "9-SNZP-m1-90_S15",
                MII_rep1    = "1-MII-B-m1-280_S1",
                MII_rep2    = "2-MII-C-m1-280_S2")

# ---------------------------- checks -----------------------------------------
for (f in c(chip_file, meth_file)) {
    if (!file.exists(f)) stop("File not found: ", f, "\nWorking directory: ", getwd(), call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(...) file.path(out_dir, paste0(...))

# ---------------------------- 1. load and normalize --------------------------
counts  <- read_bin_matrix(chip_file)
needed  <- c(fig3a_cols, k36_cols, fig3c_cols)
missing <- setdiff(needed, colnames(counts))
if (length(missing) > 0) {
    stop("Columns not found in ", chip_file, ": ", paste(missing, collapse = ", "), call. = FALSE)
}
chip <- normalize_by_top_fraction(counts[, needed], fraction = 0.05)
bins <- rownames(chip)

meth <- read_bin_matrix(meth_file)
if (!meth_col %in% colnames(meth)) meth_col <- colnames(meth)[1]
meth <- meth[match(bins, rownames(meth)), meth_col]
if (all(is.na(meth))) stop("No bins in ", meth_file, " match ", chip_file, call. = FALSE)
if (max(meth, na.rm = TRUE) > 1) meth <- meth / 100      # % -> fraction

fig3a <- setNames(chip[, fig3a_cols], names(fig3a_cols))
fig3b <- data.frame(MII_H3K36me3 = rowMeans(chip[, k36_cols]),
                    MII_DNAme    = meth, row.names = bins)
fig3c <- setNames(chip[, fig3c_cols], names(fig3c_cols))

# ---------------------------- 2. mixture modelling ---------------------------
mix_input <- data.frame(
    MII_H3K4me3  = (fig3a$MII_H3K4me3_rep1  + fig3a$MII_H3K4me3_rep2)  / 2,
    MII_H3K27me3 = (fig3a$MII_H3K27me3_rep1 + fig3a$MII_H3K27me3_rep2) / 2,
    row.names = bins
)
mix <- fit_mixture_models(mix_input, k = 2, seed = seed,
                          threshold_rule = "mean_of_means", scale = "log",
                          start = "quantile",
                          plot_file = out("mixture_model_fits.pdf"))
fit_thr <- setNames(exp(mix$table$threshold_ln), mix$table$dataset)   # normalized scale
thr <- if (is.null(fixed_thresholds)) fit_thr else fixed_thresholds[names(mix_input)]
write.csv(data.frame(dataset = names(mix_input),
                     threshold_used = unname(thr[names(mix_input)]),
                     threshold_mixture_fit = unname(fit_thr[names(mix_input)])),
          out("thresholds.csv"), row.names = FALSE)
write.csv(mix$table, out("mixture_model_parameters.csv"), row.names = FALSE)

# ---------------------------- 3. groups --------------------------------------
k4me3  <- mix_input$MII_H3K4me3  > thr["MII_H3K4me3"]
k27me3 <- mix_input$MII_H3K27me3 > thr["MII_H3K27me3"]
group  <- ifelse(k27me3 & !k4me3, "H3K27me3",
          ifelse(k27me3 &  k4me3, "Overlap",
          ifelse(!k27me3 & k4me3, "H3K4me3", "Negative")))
group_levels <- c("H3K27me3", "Overlap", "H3K4me3", "Negative")

set.seed(seed)
shuffled <- sample(seq_along(bins))
o <- shuffled[order(match(group[shuffled], group_levels))]

sizes <- as.data.frame(table(factor(group, levels = group_levels)))
colnames(sizes) <- c("group", "n_bins")
sizes$percent <- round(100 * sizes$n_bins / sum(sizes$n_bins), 1)
write.csv(sizes, out("Fig3_group_sizes.csv"), row.names = FALSE)
print(sizes, row.names = FALSE)

# ---------------------------- 4. heatmaps ------------------------------------
group_labels <- setNames(paste0(sizes$group, " (n = ", formatC(sizes$n_bins, format = "d", big.mark = ","), ")"),
                         as.character(sizes$group))
row_annot <- data.frame(group = factor(group_labels[group[o]], levels = group_labels),
                        row.names = bins[o])
annot_colors <- list(group = setNames(c("#B8A9D9", "#8C9AA8", "#9CCBAF", "white"), group_labels))

all_values <- cbind(fig3a, fig3b, fig3c)[o, ]
write.table(data.frame(bin = bins[o], group = group[o], all_values, check.names = FALSE),
            out("Fig3_bins_groups_values.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_panel <- function(values, file, title, gaps_col = NULL, width = 4) {
    hm <- as.matrix(values[o, , drop = FALSE])
    rownames(hm) <- bins[o]
    pheatmap(hm, cluster_cols = FALSE, cluster_rows = FALSE,
             color = colorRampPalette(brewer.pal(9, "Greys"))(100),
             breaks = seq(0, 1, length.out = 101),
             annotation_row = row_annot, annotation_colors = annot_colors,
             show_rownames = FALSE, gaps_col = gaps_col, angle_col = "270",
             main = title, filename = out(file), width = width, height = 10)
}
n_lab <- paste0("mm10 ", bin_size, " bins (n = ", formatC(length(bins), format = "d", big.mark = ","), ")")
plot_panel(fig3a, "Fig3a_grouped_heatmap.tiff", n_lab, gaps_col = 2)
plot_panel(fig3b, "Fig3b_K36me3_DNAme_heatmap.tiff", "MII", gaps_col = 1, width = 3.5)
plot_panel(fig3c, "Fig3c_H3K4me1_heatmap.tiff", "H3K4me1", gaps_col = c(2, 4, 6), width = 6)

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
