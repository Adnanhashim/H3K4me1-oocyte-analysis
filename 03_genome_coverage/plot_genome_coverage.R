#!/usr/bin/env Rscript
# =============================================================================
# plot_genome_coverage.R
#
# Bar plot of H3K4me1 genome coverage (%) per cell type (Fig. 1c).
# Reads genome_coverage.txt written by genome_coverage.sh.
#
# Usage (RStudio):
#   1. Open H3K4me1-oocyte-analysis.Rproj (working directory = repository root)
#   2. Copy genome_coverage.txt from the server to results/genome_coverage/
#   3. Check the paths under "input / output paths" below, then click "Source"
#
# Usage (command line, optional):
#   Rscript plot_genome_coverage.R <genome_coverage.txt> <labels.tsv> <output_prefix>
#
#   genome_coverage.txt  output of genome_coverage.sh
#   labels.tsv           columns 'sample' and 'label'; sets the bar labels and
#                        their order (see genome_coverage_labels.tsv)
#   output_prefix        e.g. figures/Fig1/Fig1c_genome_coverage
#
# Output:
#   <output_prefix>.svg                 the plot (6 x 4 in)
#   <output_prefix>_plotted_values.csv  the values shown in the plot
#
# Requires: R >= 4.0, ggplot2, svglite
# =============================================================================

# Load required library
suppressPackageStartupMessages(library(ggplot2))

# ---------------------------- input / output paths ---------------------------
# Relative to the repository root (RStudio: open H3K4me1-oocyte-analysis.Rproj).
coverage_file <- "results/genome_coverage/genome_coverage.txt"
labels_file   <- "03_genome_coverage/genome_coverage_labels.tsv"
out_prefix    <- "figures/Fig1/Fig1c_genome_coverage"

# Command-line arguments override the paths above (Rscript use only)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3) {
    coverage_file <- args[1]
    labels_file   <- args[2]
    out_prefix    <- args[3]
} else if (length(args) > 0) {
    stop("Usage: Rscript plot_genome_coverage.R <genome_coverage.txt> <labels.tsv> <output_prefix>",
         call. = FALSE)
}

for (f in c(coverage_file, labels_file)) {
    if (!file.exists(f)) {
        stop("File not found: ", f, "\nWorking directory is: ", getwd(),
             "\n(Open H3K4me1-oocyte-analysis.Rproj so paths are relative to the repository root.)",
             call. = FALSE)
    }
}

# ---------------------------- define data ------------------------------------
coverage <- read.delim(coverage_file, stringsAsFactors = FALSE)
labels   <- read.delim(labels_file, stringsAsFactors = FALSE, comment.char = "#")

if (!all(c("sample", "label") %in% colnames(labels))) {
    stop(labels_file, " must be TAB-separated with columns 'sample' and 'label' ",
         "(spaces instead of tabs is the usual cause).", call. = FALSE)
}

missing <- setdiff(labels$sample, coverage$sample)
if (length(missing) > 0) {
    stop("Samples in ", labels_file, " not found in ", coverage_file, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
}

# Create data frame (bar order = order in labels file)
df <- data.frame(
    cell_types      = factor(labels$label, levels = labels$label),
    genome_coverage = coverage$genome_coverage_pct[match(labels$sample, coverage$sample)]
)

# ---------------------------- create plot ------------------------------------
p <- ggplot(data = df, aes(x = cell_types, y = genome_coverage)) +
    geom_bar(stat = "identity", fill = "blue", color = "black", width = 0.5) +
    geom_text(aes(label = genome_coverage), vjust = -0.3, size = 3.5) +
    labs(x = "Cell Type",
         y = "Genome Coverage (%)",
         title = "Genome Coverage Percentage by Cell Type") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )

# ---------------------------- save plot --------------------------------------
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

ggsave(paste0(out_prefix, ".svg"), plot = p, width = 6, height = 4, units = "in", dpi = 300)
write.csv(df, paste0(out_prefix, "_plotted_values.csv"), row.names = FALSE)

message("Saved: ", out_prefix, ".svg and ", out_prefix, "_plotted_values.csv")
print(df, row.names = FALSE)
