# =============================================================================
# Fig4d_human_grouped_heatmap.R
#
# Fig. 4d  Grouped heatmap of human oocyte H3K4me1 (GO1, GO2, MI/II), 4-cell
#          H3K4me3 (rep1, rep2), GV H3K27me3 and MII DNA methylation in hg38
#          10 kb bins, grouped by H3K27me3 (GV) and H3K4me3 (4-cell) status.
#
# Method
#   1. Read counts in hg38 10 kb bins, each dataset normalized to its top 5%
#      signal (00_shared/normalization.R). GO1 / GO2 H3K4me1 = mean of two
#      oocytes each; 4-cell H3K4me3 = mean of two replicates for grouping.
#   2. Positive/negative thresholds by mixture modelling
#      (05_mixture_model/mixture_model.R): two-component normal mixture on the
#      normalized signal; threshold = mean of the two component means.
#      The figure uses the thresholds of the original fit (mixtools
#      normalmixEM, random starting values): MI/II H3K4me1 0.60, 4-cell
#      H3K4me3 0.3946, GV H3K27me3 0.408. The script refits the models with a
#      fixed seed and reports both in thresholds.csv.
#   3. Bins grouped as: H3K27me3+ H3K4me3-  |  H3K27me3+ H3K4me3+ (overlap)  |
#      H3K27me3- H3K4me3+  |  negative. Random (seeded) order within groups.
#   4. Heatmap (0-1, white-black); MII DNA methylation shown as fraction; the
#      binary mixture-model calls (4-cell H3K4me3, GV H3K27me3) as the last
#      two columns.
#   The H3K4me3+/H3K27me3- bins are used for the subclustering in Fig. 4e
#   (06_kmeans/Fig4e_human_subclustering_400bp.R).
#
# Input (data/Fig4/)
#   human_10kb_GV_H3K27me3_4C_H3K4me3_counts.txt     10 kb read counts
#   human_10kb_MI-MII_H3K4me1_counts.txt                       10 kb read counts
#   human_10kb_GO1_GO2_H3K4me1_counts.txt    10 kb read counts (GO1/GO2)
#   human_10kb_MII_5mC.txt                            10 kb DNA methylation (%)
# Output (figures/Fig4/Fig4d/)
#   Fig4d_grouped_heatmap.tiff        Fig. 4d
#   Fig4d_bins_groups_values.txt      bin, group, plotted values and binary calls
#                                     (heatmap order)
#   Fig4d_group_sizes.csv             number of bins per group
#   thresholds.csv                    thresholds used (original fit) and refitted
#   mixture_model_parameters.csv      refitted mixture-model parameters
#   mixture_model_fits.pdf            fitted mixture per dataset
#   K4me3pos_K27me3neg_bins.txt       10 kb bins used for Fig. 4e
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
if (!"start_sd" %in% names(formals(fit_mixture_models))) {
    stop("05_mixture_model/mixture_model.R is an older version - replace it with the current one.",
         call. = FALSE)
}

# ---------------------------- input / output ---------------------------------
data_dir   <- "data/Fig4"
gv_4c_file <- file.path(data_dir, "human_10kb_GV_H3K27me3_4C_H3K4me3_counts.txt")
mi_file    <- file.path(data_dir, "human_10kb_MI-MII_H3K4me1_counts.txt")
go_file    <- file.path(data_dir, "human_10kb_GO1_GO2_H3K4me1_counts.txt")
meth_file  <- file.path(data_dir, "human_10kb_MII_5mC.txt")
out_dir    <- "figures/Fig4/Fig4d"
seed       <- 123

# Thresholds (normalized signal) from the original mixture-model fit
original_thresholds <- c(MI_II_H3K4me1 = 0.60, C4_H3K4me3_avg = 0.3946, GV_H3K27me3 = 0.408)

# ---------------------------- checks -----------------------------------------
for (f in c(gv_4c_file, mi_file, go_file, meth_file)) {
    if (!file.exists(f)) stop("File not found: ", f, "\nWorking directory: ", getwd(), call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(...) file.path(out_dir, paste0(...))

# ---------------------------- 1. load and normalize --------------------------
gv_4c <- read_bin_matrix(gv_4c_file)
mi    <- read_bin_matrix(mi_file)
bins  <- union(rownames(gv_4c), rownames(mi))
chip  <- cbind(gv_4c[bins, , drop = FALSE], mi[bins, , drop = FALSE])
rownames(chip) <- bins
chip  <- normalize_by_top_fraction(chip, fraction = 0.05)

go <- normalize_by_top_fraction(read_bin_matrix(go_file), fraction = 0.05)
go <- go[match(bins, rownames(go)), ]

meth <- read_bin_matrix(meth_file)
meth <- meth[match(bins, rownames(meth)), 1] / 100       # % -> fraction

values <- data.frame(
    GO1_H3K4me1   = (go[["13-hPriA-K4m1"]] + go[["14-hPriB-K4m1"]]) / 2,
    GO2_H3K4me1   = (go[["15-hSecA-K4m1"]] + go[["16-hSecB-K4m1"]]) / 2,
    MI_II_H3K4me1 = chip[["K4me1_MI-II"]],
    C4_H3K4me3_rep1 = chip[["4cell_2PN_H3K4me3_rep1"]],
    C4_H3K4me3_rep2 = chip[["4cell_2PN_H3K4me3_rep2"]],
    GV_H3K27me3   = chip[["GV_H3K27me3"]],
    MII_DNAme     = meth,
    row.names = bins
)

# ---------------------------- 2. mixture modelling ---------------------------
mix_input <- data.frame(
    MI_II_H3K4me1  = values$MI_II_H3K4me1,
    C4_H3K4me3_avg = (values$C4_H3K4me3_rep1 + values$C4_H3K4me3_rep2) / 2,
    GV_H3K27me3    = values$GV_H3K27me3,
    row.names = bins
)
mix <- fit_mixture_models(mix_input, k = 2, seed = seed,
                          threshold_rule = "mean_of_means", scale = "linear",
                          start = "quantile", start_sd = c(0.05, 0.1, 0.2, 0.5),
                          plot_file = out("mixture_model_fits.pdf"))
refit <- setNames(exp(mix$table$threshold_ln), mix$table$dataset)   # normalized scale
write.csv(data.frame(dataset = names(original_thresholds),
                     threshold_used = original_thresholds,
                     threshold_refit = refit[names(original_thresholds)]),
          out("thresholds.csv"), row.names = FALSE)
write.csv(mix$table, out("mixture_model_parameters.csv"), row.names = FALSE)
thr <- original_thresholds

# ---------------------------- 3. groups --------------------------------------
k4me3  <- mix_input$C4_H3K4me3_avg > thr["C4_H3K4me3_avg"]
k27me3 <- mix_input$GV_H3K27me3    > thr["GV_H3K27me3"]
group  <- ifelse(k27me3 & !k4me3, "H3K27me3",
          ifelse(k27me3 &  k4me3, "Overlap",
          ifelse(!k27me3 & k4me3, "H3K4me3", "Negative")))
group_levels <- c("H3K27me3", "Overlap", "H3K4me3", "Negative")

set.seed(seed)
shuffled <- sample(seq_along(bins))
o <- shuffled[order(match(group[shuffled], group_levels))]

sizes <- as.data.frame(table(factor(group, levels = group_levels)))
colnames(sizes) <- c("group", "n_bins")
write.csv(sizes, out("Fig4d_group_sizes.csv"), row.names = FALSE)
print(sizes, row.names = FALSE)

writeLines(bins[group == "H3K4me3"], out("K4me3pos_K27me3neg_bins.txt"))

# ---------------------------- 4. heatmap -------------------------------------
# Group labels with the number of bins, e.g. "H3K4me3 (n = 56,589)"
group_labels <- setNames(paste0(sizes$group, " (n = ", formatC(sizes$n_bins, format = "d", big.mark = ","), ")"),
                         as.character(sizes$group))
# Binary (mixture-model) calls shown as the last two columns: 1 = positive
values$Binary_C4_H3K4me3 <- as.integer(k4me3)
values$Binary_GV_H3K27me3 <- as.integer(k27me3)
hm <- as.matrix(values[o, ])
write.table(data.frame(bin = bins[o], group = group[o], hm, check.names = FALSE),
            out("Fig4d_bins_groups_values.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

pheatmap(hm, cluster_cols = FALSE, cluster_rows = FALSE,
         color = colorRampPalette(brewer.pal(9, "Greys"))(100),
         breaks = seq(0, 1, length.out = 101),
         annotation_row = data.frame(group = factor(group_labels[group[o]], levels = group_labels),
                                     row.names = rownames(hm)),
         show_rownames = FALSE, gaps_col = c(3, 5, 6, 7),
         main = "hg38 10 kb bins, top 5% normalized",
         filename = out("Fig4d_grouped_heatmap.tiff"), width = 6, height = 10)

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
