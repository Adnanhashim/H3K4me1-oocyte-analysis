# =============================================================================
# Fig1e_1f_K4me1_vs_K27ac.R
#
# Fig. 1e  2D histograms of H3K4me1 vs H3K27ac signal in TSS-distal 10 kb bins
#          for P10 growing oocytes (GO) and ENCODE cell types/tissues.
# Fig. 1f  Percentage of TSS-distal H3K4me1-positive 10 kb bins that are
#          H3K27ac-negative.
#
# Method
#   1. Read counts in 10 kb bins (mm10); P10 GO H3K4me1 = mean of the two
#      biological replicates.
#   2. Each dataset is normalized to its 95th percentile (all bins).
#   3. TSS-distal bins: 10 kb bins whose nearest TSS (refTSS v3.3) is more
#      than 2 kb away.
#   4. Fig. 1e: 2D histogram of ln-transformed signal; Pearson r on the
#      normalized signal.
#   5. Positive/negative thresholds per dataset by mixture modelling
#      (05_mixture_model/mixture_model.R): two-component normal mixture fitted
#      to the ln-transformed signal of bins with signal; threshold = mean of
#      the two component means. Shown as dashed lines in Fig. 1e.
#   6. Fig. 1f: % of H3K4me1-positive bins that are H3K27ac-negative. Bins
#      without signal are negative.
#
# Input (data/Fig1/)
#   Fig1F_mm10_10KB_bins_scatterplot.txt            all 10 kb bins
#   Fig1F_mm10_10KB_bins_TSS_distalpeaks_input.txt  TSS-distal 10 kb bins
#
# Output (figures/Fig1/Fig1e_1f/)
#   Fig1e_2Dhistograms.pdf/.svg       Fig. 1e
#   Fig1f_proportions.pdf/.svg        Fig. 1f
#   Fig1f_proportions.csv             bin counts and % per cell type
#   pearson_r.csv                     r per cell type
#   thresholds.csv                    mixture-model parameters and thresholds
#   mixture_model_fits.pdf            fitted mixture per dataset
#   distal_bins_normalized.txt        normalized signal of the plotted bins
#   sessionInfo.txt                   R and package versions
#
# The final figure layout was edited in Adobe Illustrator.
#
# Usage (RStudio): open H3K4me1-oocyte-analysis.Rproj and click "Source".
# Requires: ggplot2, hexbin, gridExtra, svglite, mixtools, fBasics
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(gridExtra)
    library(svglite)
})
source("05_mixture_model/mixture_model.R")

# ---------------------------- input / output ---------------------------------
data_dir         <- "data/Fig1"
all_bins_file    <- file.path(data_dir, "Fig1F_mm10_10KB_bins_scatterplot.txt")
distal_bins_file <- file.path(data_dir, "Fig1F_mm10_10KB_bins_TSS_distalpeaks_input.txt")
out_dir          <- "figures/Fig1/Fig1e_1f"

# ---------------------------- parameters -------------------------------------
cell_types  <- c("p10", "ES", "Kidney", "Liver", "Smint", "Testis")
cell_labels <- c(p10 = "P10 GO", ES = "ESC", Kidney = "Kidney", Liver = "Liver",
                 Smint = "Small intestine", Testis = "Testis")
seed <- 123

# ---------------------------- checks -----------------------------------------
for (f in c(all_bins_file, distal_bins_file)) {
    if (!file.exists(f)) {
        stop("File not found: ", f, "\nWorking directory: ", getwd(),
             "\n(Open H3K4me1-oocyte-analysis.Rproj so paths are relative to the repository root.)",
             call. = FALSE)
    }
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(...) file.path(out_dir, paste0(...))

# ---------------------------- 1-2. load and normalize ------------------------
df <- read.table(all_bins_file, header = TRUE, sep = "\t")
df$Bins <- paste(df$chr, df$sp, df$ep, sep = "_")
df$p10_K4me1 <- (df$p10A_K4me1 + df$p10B_K4me1) / 2

value_cols <- as.vector(t(outer(cell_types, c("K4me1", "K27ac"), paste, sep = "_")))
norm_all <- cbind(df["Bins"], normalize_by_quantile(df[, value_cols], probs = 0.95))

# ---------------------------- 3. TSS-distal bins -----------------------------
distal <- read.table(distal_bins_file, header = TRUE, sep = "\t")
distal_bins <- unique(paste(distal$chr, distal$sp, distal$ep, sep = "_"))

bins <- norm_all[norm_all$Bins %in% distal_bins, ]
rownames(bins) <- bins$Bins
bins <- bins[, value_cols]
write.table(bins, out("distal_bins_normalized.txt"), sep = "\t", quote = FALSE)
message("Bins: ", nrow(norm_all), " total, ", nrow(bins), " TSS-distal")

# ---------------------------- 5. mixture modelling ---------------------------
mix <- fit_mixture_models(bins, k = 2, seed = seed,
                          threshold_rule = "mean_of_means", scale = "log", start = "quantile",
                          plot_file = out("mixture_model_fits.pdf"))
if (!all(c("dataset", "threshold_ln") %in% colnames(mix$table))) {
    stop("mixture_model.R is an older version - replace it with the current one from the repository.",
         call. = FALSE)
}
write.csv(mix$table, out("thresholds.csv"), row.names = FALSE)
thresholds <- setNames(mix$table$threshold_ln, mix$table$dataset)

# ---------------------------- 4. Fig. 1e -------------------------------------
pearson_r <- sapply(cell_types, function(ct) {
    cor(bins[[paste0(ct, "_K4me1")]], bins[[paste0(ct, "_K27ac")]])
})
write.csv(data.frame(cell_type = cell_types, pearson_r = pearson_r),
          out("pearson_r.csv"), row.names = FALSE)

theme_fig <- theme(
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
)

plot_2d <- function(ct) {
    d <- data.frame(x = log(bins[[paste0(ct, "_K4me1")]]),
                    y = log(bins[[paste0(ct, "_K27ac")]]))
    d <- d[is.finite(d$x) & is.finite(d$y), ]
    ggplot(d, aes(x, y)) +
        geom_hex(bins = 120, aes(fill = after_stat(count))) +
        scale_fill_gradientn(name = "Count",
                             colors = colorRampPalette(c("lightblue", "darkblue"))(9),
                             limits = c(1, 500)) +
        geom_vline(xintercept = thresholds[paste0(ct, "_K4me1")], linetype = "dashed") +
        geom_hline(yintercept = thresholds[paste0(ct, "_K27ac")], linetype = "dashed") +
        annotate("text", x = Inf, y = Inf, label = paste("r =", round(pearson_r[ct], 2)),
                 hjust = 1, vjust = 1, size = 4) +
        labs(title = cell_labels[ct], x = "ln H3K4me1", y = "ln H3K27ac") +
        theme_fig
}

fig1e <- arrangeGrob(grobs = lapply(cell_types, plot_2d), ncol = 2, nrow = 3)
for (ext in c("pdf", "svg")) {
    ggsave(out("Fig1e_2Dhistograms.", ext), fig1e, width = 14, height = 21)
}

# ---------------------------- 6. Fig. 1f -------------------------------------
bins_ln <- log_transform(bins, zero_value = -Inf)

fig1f_table <- do.call(rbind, lapply(cell_types, function(ct) {
    k4 <- bins_ln[[paste0(ct, "_K4me1")]] >= thresholds[paste0(ct, "_K4me1")]
    ac <- bins_ln[[paste0(ct, "_K27ac")]] >= thresholds[paste0(ct, "_K27ac")]
    data.frame(cell_type = ct,
               K4me1pos_K27acneg = sum(k4 & !ac),
               K4me1pos_K27acpos = sum(k4 & ac),
               K4me1neg_K27acpos = sum(!k4 & ac),
               K4me1neg_K27acneg = sum(!k4 & !ac))
}))
fig1f_table$percent <- 100 * fig1f_table$K4me1pos_K27acneg /
    (fig1f_table$K4me1pos_K27acneg + fig1f_table$K4me1pos_K27acpos)
write.csv(fig1f_table, out("Fig1f_proportions.csv"), row.names = FALSE)

fig1f <- ggplot(fig1f_table, aes(x = factor(cell_labels[cell_type], levels = cell_labels),
                                 y = percent)) +
    geom_bar(stat = "identity", fill = "grey50", width = 0.5) +
    labs(title = "H3K4me1-positive bins that are H3K27ac-negative",
         x = NULL, y = "Proportion (%)") +
    ylim(0, 100) +
    theme_classic()
for (ext in c("pdf", "svg")) {
    ggsave(out("Fig1f_proportions.", ext), fig1f, width = 8, height = 6)
}

# ---------------------------- summary ----------------------------------------
print(data.frame(cell_type = cell_types, r = round(pearson_r, 2),
                 percent_K27ac_neg = round(fig1f_table$percent, 1)), row.names = FALSE)
writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
