#' Fit cosme model
#' @description
#' A short description...
#'
# TODO better description of parameters
#' @param model a `cosme_model` object
#' @param data a `data.frame`
#' @param option a `list` object with options passed to various methods. See details
#' @details
#' Additional details...
#'
#' @import cli
#' @export
cosme <- function(model, data, option) {
    # Check that it is a model
    if (!is.model(model)) model <- make_model(model)

    # Check for options
    if (missing(option)) option <- list()

    # Create input list
    inp <- list(model = model, data = data)

    # Run Frequentist
    cli::cli_alert('Fitting Frequentist')
    flush.console()
    freq  <- do.call(fit_freq,  c(inp, option$freq))

    # Run Bayesian
    cli::cli_alert('Fitting Bayesian')
    flush.console()
    cli_progress_bar(
        total = 2,
        clear = TRUE,
        format = paste0(
            "    {bmodel} {cli::pb_current}/{cli::pb_total} ",
            "| Elapsed: {cli::pb_elapsed}"
        )
    )
    bmodel <- 'Full:'
    cli_progress_update(set = 1)
    bout <- capture.output(
        bayes <- do.call(fit_bayes, c(inp, option$bayes)),
        nullfile()
    )
    bmodel <- 'Baseline:'
    cli_progress_update(set = 2)
    bbout <- capture.output(
        baseline <- do.call(
            fit_bayes, list(
                model = baseline_model(freq), data = data, option$bayes
            )
        ),
        nullfile()
    )
    flush.console()
    cli_progress_done()

    # Run information theory
    cli::cli_alert('Fitting Information Theory')
    flush.console()
    info <- suppressWarnings(do.call(fit_info, c(inp, option$info)))

    # Create `cosme_fit` and return
    structure(
        list2env(
            list(
                freq  = freq,
                bayes = bayes,
                info  = info
            ),
            parent = emptyenv()
        ),
        class = c("cosme_fit", "environment"),
        b_base = baseline
    )
}

#' Check if `cosme_fit`
#' @noRd
is.fit <- function(x) {
    inherits(x, "cosme_fit")
}

#' Obtain Estimates from `comse_fit`
#' @description
#' A short description...
#'
# TODO better description of parameters
#' @param object a `cosme_model` object
#' @import lavaan blavaan
#' @export
estimates <- function(object) {
    # Check that it is a `cosme_fit`
    if (!is.fit(object)) {
        throw_error(c(
            "x" = "The object is not of class {.cls cosme_fit}"
        ))
    }

    ## Get parameter estimates
    f <- parameterEstimates(object$freq)
    b <- parameterEstimates(object$bayes)
    i <- vapply(object$info$mod_list[[1]], \(.) .@ParTable$est, numeric(NROW(f)))

    ## Extract Bayesian posterior intervals
    # Have to capture output to prevent printing
    tmp <- capture.output(x <- summary(object$bayes), nullfile())
    # Extract intervals |> make numbers |> make matrix -> ci
    x[, c('pi.lower', 'pi.upper')] |> as.numeric() |> matrix(ncol = 2) -> bpi
    # Change NA to estimate
    is.na(bpi) |> ifelse(b$est, bpi) -> bpi
    # blavaan::blavInspect(o$bayes, 'hpd') # Doesnt match same rows columns
    ## Extract information theory interval
    iint <- t(apply(i, 1, quantile, prob = c(0.025, 0.975)))

    ## Return
    structure(
        list(
            estimate = cbind(
                freq  = f$est,
                bayes = b$est,
                info  = apply(i, 1, median)
            ),
            interval_lower = cbind(
                freq  = f$ci.lower,
                bayes = bpi[,1],
                info  = iint[,1]
            ),
            interval_upper = cbind(
                freq  = f$ci.upper,
                bayes = bpi[,2],
                info  = iint[,2]
            )
        ),
        pname = apply(
            f[,c('lhs', 'op', 'rhs')],
            1, paste0, collapse = ' '
        )
    )
}


#' Obtain fit statistics from `comse_fit`
#' @description
#' A short description...
#'
# TODO better description of parameters
#' @param object a `cosme_model` object
#' @import semTools lavaan blavaan
#' @export
fit <- function(object) {
    # Check that it is a `cosme_fit`
    if (!is.fit(object)) {
        throw_error(c(
            "x" = "The object is not of class {.cls cosme_fit}"
        ))
    }

    # Obtain bayesian fit
    # Takes awhile
    bfit <- blavaan::blavFitIndices(
        o$bayes, baseline.model = attr(o, 'b_base'),
        fit.measures = c('BRMSEA', 'BCFI')
    )


    # Return output
    list(
        # DOES NOT CURRENTLY WORK DUE TO BUG IN semTools
        # Need to wait for semTools to update
        # freq  = list(lavaan::fitmeasures(object$freq), semTools::moreFitIndices(object$freq)),
        freq  = lavaan::fitmeasures(object$freq, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr")),
        bayes =c(
            ppp = lavaan::fitMeasures(o$bayes, c("ppp"))[[1]],
            rmsea = mean(bfit@indices$BRMSEA),
            cfi = mean(bfit@indices$BCFI)
        ),
        info  = o$info$fit_list[[1]]
    )
}
