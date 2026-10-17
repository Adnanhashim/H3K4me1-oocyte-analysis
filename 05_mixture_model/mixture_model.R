# =============================================================================
# mixture_model.R
#
# Shared functions for the Gaussian mixture modelling used to classify 10 kb
# bins as positive or negative for a histone mark (approach adapted from
# Dileep & Gilbert). All 10 kb-bin mixture modelling in the manuscript uses
# these functions; each figure script sources this file and sets its own
# parameters:
#
#   Figure   Bins                    Normalization            k  Threshold
#   1e/1f    10 kb, TSS-distal       / 95th percentile        2  mean of the two
#                                    (per dataset)               component means
#                                                                (fit on ln values)
#   4b       10 kb, bosTau8          top 5% (per dataset)     3  mean of the two
#                                                                lowest component means
#                                                                (fit on normalized values)
#   4d       10 kb, hg38             top 5% (per dataset)     2  mean of the two
#                                                                component means
#                                                                (fit on normalized values)
#
# Per dataset (column):
#   1. signal normalized (00_shared/normalization.R)
#   2. mixtools::normalmixEM fitted with k components
#   3. threshold = mean (or median) of the fitted component means
#   4. density plot of the data with the fitted components
#
# Usage: source("05_mixture_model/mixture_model.R") from a figure script.
# Requires: mixtools, fBasics
# =============================================================================

suppressPackageStartupMessages({
    library(mixtools)
    library(fBasics)
})

source("00_shared/normalization.R")   # normalize_by_quantile, normalize_by_top_fraction, log_transform

# ---------------------------- plotting ---------------------------------------

#' Density of the data with the fitted components and the component means.
plot_mixture <- function(x, model, title, density_adjust = 0.5, xlim = NULL) {
    d  <- density(x, adjust = density_adjust)
    xa <- min(d$x) - 0.1; xb <- max(d$x) + 0.1
    if (!is.null(xlim)) { xa <- xlim[1]; xb <- xlim[2] }
    yb <- max(d$y) + 0.1

    plot(d, xlim = c(xa, xb), ylim = c(0, yb), main = title)
    abline(v = model$mu)

    comp_cols <- c("red", "green", "blue", "orange", "purple")
    xs <- seq(xa, xb, length.out = 1000)
    for (i in seq_along(model$mu)) {
        lines(xs, model$lambda[i] * dnorm(xs, model$mu[i], model$sigma[i]), col = comp_cols[i])
    }
}

# ---------------------------- fitting ----------------------------------------

#' Fit a k-component normal mixture to every column of 'mat'.
#'
#' @param mat             data.frame of normalized signal (bins x datasets)
#' @param k               number of components
#' @param maxit           maximum EM iterations
#' @param seed            random seed
#' @param threshold_rule  "mean_of_means": mean of all component means (k = 2)
#'                        "low_mid_mean" : mean of the two lowest component
#'                                         means (k = 3; the third component
#'                                         takes the few very high bins)
#'                        "middle_mean"  : median of the component means
#' @param scale           "log": fit on ln values of bins with signal
#'                        "linear": fit on the normalized values
#' @param start           "quantile": fixed starting values (means at fixed
#'                        quantiles, start_sd as SDs), the fit
#'                        with the highest log-likelihood is kept
#'                        "random": normalmixEM random starting values
#' @param plot_file       multi-page PDF with one fit per page (NULL = none)
#' @param density_adjust  bandwidth adjustment for the density plots
#' @param start_sd        starting SDs tried with start = "quantile"
#' @param plot_xlim       x-axis range of the fit plots (NULL = data range)
#' @return list(table = one row per column incl. threshold on ln scale,
#'              models = fitted normalmixEM objects)
fit_mixture_models <- function(mat, k = 2, maxit = 100000, seed = 123,
                               threshold_rule = c("mean_of_means", "low_mid_mean", "middle_mean"),
                               scale = c("log", "linear"),
                               start = c("quantile", "random"),
                               plot_file = NULL, density_adjust = 0.5,
                               start_sd = c(0.3, 0.7, 1.2), plot_xlim = NULL) {
    threshold_rule <- match.arg(threshold_rule)
    scale <- match.arg(scale)
    start <- match.arg(start)
    if (!is.null(plot_file)) {
        pdf(plot_file, width = 6, height = 8, paper = "special")
        on.exit(dev.off(), add = TRUE)
    }

    models <- list()
    rows   <- list()

    for (columnname in colnames(mat)) {
        message("Mixture model (k = ", k, "): ", columnname)
        x <- mat[[columnname]]
        if (scale == "log") { x <- log(x); x <- x[is.finite(x)] }

        set.seed(seed)
        model <- NULL
        if (start == "quantile") {
            probs <- if (k == 2) c(0.25, 0.75) else if (k == 3) c(0.25, 0.75, 0.999) else
                     seq(0.1, 0.9, length.out = k)
            for (s0 in start_sd) {
                m <- tryCatch(
                    suppressMessages(normalmixEM(x = x, k = k, maxit = maxit,
                                                 mu = unname(quantile(x, probs)),
                                                 sigma = rep(s0, k), lambda = rep(1 / k, k))),
                    error = function(e) NULL)
                if (!is.null(m) && (is.null(model) || m$loglik > model$loglik)) model <- m
            }
        } else {
            model <- tryCatch(suppressMessages(normalmixEM(x = x, k = k, maxit = maxit)),
                              error = function(e) NULL)
        }
        if (is.null(model)) {
            warning("normalmixEM did not converge for ", columnname, "; skipped.", call. = FALSE)
            next
        }
        models[[columnname]] <- model

        threshold <- switch(threshold_rule,
                            mean_of_means = mean(model$mu),
                            low_mid_mean  = mean(sort(model$mu)[1:2]),
                            middle_mean   = median(model$mu))

        rows[[columnname]] <- data.frame(
            dataset          = columnname,
            skewness         = skewness(x)[1],
            low_mean         = min(model$mu),
            high_mean        = max(model$mu),
            threshold_ln     = if (scale == "log") threshold else log(threshold),
            scale            = scale,
            loglik           = model$loglik,
            mu               = paste(signif(model$mu, 6), collapse = ";"),
            sigma            = paste(signif(model$sigma, 6), collapse = ";"),
            lambda           = paste(signif(model$lambda, 6), collapse = ";"),
            stringsAsFactors = FALSE
        )

        if (!is.null(plot_file)) plot_mixture(x, model, columnname, density_adjust, plot_xlim)
    }

    tab <- do.call(rbind, rows)
    rownames(tab) <- NULL
    list(table = tab, models = models)
}
