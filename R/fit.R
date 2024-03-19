#' Fit Frequentist model via `lavaan`
#' @description
#' A short description...
#'
#' @param model a `cosme_model` object
#' @param data a `data.frame` object
#' @param ... passed to `lavaan::cfa`
#' @importFrom lavaan cfa lavaan
#' @noRd
fit_freq <- function(model, data, ...) {
    model |> as.lavaan() |> cfa(data, ...)
}

#' Fit Bayesian model via `blavaan`
#' @description
#' A short description...
#'
#' @param model a `cosme_model` object
#' @param data a `data.frame` object
#' @param ... passed to `blavaan::bcfa`
#' @importFrom blavaan bcfa blavaan
#' @noRd
fit_bayes <- function(model, data, ...) {
    model |> as.lavaan() |> bcfa(data, ...)
}

#' Fit Information Theory model via implementation of `ockhamSEM`
#' @description
#' A short description...
#'
#' @param model a `cosme_model` object
#' @param data a `data.frame` object
#' @param reps number of replicated data sets
#' @param ... passed to `lavaan::cfa`
# TODO better description of ...
#' @noRd
fit_info <- function(model, data, reps, ...) {
    # TODO print warning with missing data
    # TODO deal with reps
    if (missing(reps)) reps <- 25
    tmpmat <- diag(NCOL(data))
    colnames(tmpmat) <- rownames(tmpmat) <- names(data)

    (model
        |> as.lavaan()
        |> cfa(sample.cov = tmpmat, sample.nobs = NROW(data))
        |> run_fitprop(
            fit.measure =  c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
            rmethod = "onion",
            reps = reps,
            ...
        )
    )
}
