# =============================================================================
# genomic_element_distribution.R
#
# Distribution of broad peaks over genomic elements (promoter, gene body,
# downstream, intergenic, ...) with ChIPpeakAnno::genomicElementDistribution.
# Used for Fig. 1d (P10 GO H3K4me1 vs ENCODE cell types/tissues).
#
# Annotation : TxDb.Mmusculus.UCSC.mm10.knownGene
# Promoter   : 2000 bp upstream to 500 bp downstream of the TSS
# Downstream : 0 to 2000 bp downstream of the gene end
#
# Input : all *.broadPeak files in 'peaks_dir'
# Output (in 'out_dir'), per peak file:
#   <sample>.eps                   genomic element distribution plot
#   <sample>_plotted_values.csv    values shown in the plot (when available)
# plus sessionInfo.txt (R and package versions, incl. TxDb annotation version)
#
# Usage (RStudio):
#   1. Open H3K4me1-oocyte-analysis.Rproj (working directory = repository root)
#   2. Check the paths under "input / output paths" below, then click "Source"
#
# Requires (Bioconductor): GenomicRanges, ChIPpeakAnno,
#                          TxDb.Mmusculus.UCSC.mm10.knownGene
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("GenomicRanges", "ChIPpeakAnno", "TxDb.Mmusculus.UCSC.mm10.knownGene"))
# =============================================================================

# Load necessary libraries
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(ChIPpeakAnno)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

# ---------------------------- input / output paths ---------------------------
# Relative to the repository root (RStudio: open H3K4me1-oocyte-analysis.Rproj).
peaks_dir <- "peaks"
out_dir   <- "figures/Fig1/Fig1d_genomic_element_distribution"

# ---------------------------- parameters (as used in the manuscript) ---------
txdb            <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter_region <- c(upstream = 2000, downstream = 500)
gene_downstream <- c(upstream = 0, downstream = 2000)

# ---------------------------- input files ------------------------------------
if (!dir.exists(peaks_dir)) {
    stop("Peak directory not found: ", peaks_dir, "\nWorking directory is: ", getwd(),
         "\n(Open H3K4me1-oocyte-analysis.Rproj so paths are relative to the repository root.)",
         call. = FALSE)
}

# Get the broadPeak filenames in the directory
filenames <- list.files(peaks_dir, pattern = "\\.broadPeak$", full.names = TRUE)
if (length(filenames) == 0) stop("No .broadPeak files found in ", peaks_dir, call. = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Found ", length(filenames), " broadPeak file(s) in ", peaks_dir)

# ---------------------------- process each file ------------------------------
for (filename in filenames) {
    sample <- sub("\\.broadPeak$", "", basename(filename))
    message("Processing: ", sample)

    # Load the file into a GRanges object
    peaks <- toGRanges(filename, format = "broadPeak", header = FALSE)

    # Generate the plot and save to an EPS file
    outfilename <- file.path(out_dir, paste0(sample, ".eps"))
    setEPS()
    postscript(outfilename)
    res <- genomicElementDistribution(peaks,
                                      TxDb = txdb,
                                      promoterRegion = promoter_region,
                                      geneDownstream = gene_downstream)
    dev.off()

    # Export the plotted values (from the returned ggplot object, if present)
    plot_obj <- if (is.list(res) && !is.null(res$plot)) res$plot else NULL
    if (inherits(plot_obj, "ggplot")) {
        values <- tryCatch(ggplot2::layer_data(plot_obj), error = function(e) NULL)
        if (!is.null(values)) {
            values <- values[, vapply(values, is.atomic, logical(1)), drop = FALSE]
            write.csv(values, file.path(out_dir, paste0(sample, "_plotted_values.csv")),
                      row.names = FALSE)
        }
    }

    message("  saved: ", outfilename)
}

# ---------------------------- record versions --------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
