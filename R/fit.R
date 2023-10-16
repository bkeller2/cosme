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
    model |>
        as.lavaan() |>
        cfa(data, ...)
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
    model |>
        as.lavaan() |>
        bcfa(data, ...)
}

#' Fit Information Theory model via implementation of `ockhamSEM`
#' @description
#' A short description...
#'
#' @param model a `cosme_model` object
#' @param data a `data.frame` object
#' @param ... passed to `lavaan::cfa`
# TODO better description of ...
#' @noRd
fit_info <- function(model, data, ...) {
    NULL # TODO add this
}
