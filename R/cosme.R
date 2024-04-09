#' Fit cosme model
#' @description
#' A short description...
#'
# TODO better description of parameters
#' @param model a `cosme_model` object
#' @param data a `data.frame`
#' @param info a [base::logical] If `TRUE` then fit information theory
#' @param option a `list` object with options passed to various methods. See details
#' @details
#' Additional details...
#'
#' @import cli
#' @export
cosme <- function(model, data, info = FALSE, option) {
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
                model = baseline_model(freq),
                data = data, option$bayes
            )
        ),
        nullfile()
    )
    flush.console()
    cli_progress_done()

    # Run information theory
    if (info) {
        cli::cli_alert('Fitting Information Theory')
        flush.console()
        info <- suppressWarnings(do.call(fit_info, c(inp, option$info)))
    } else {
        info <- NULL
    }

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

#' Summarize from `comse_fit`
#' @description
#' A short description...
#'
# TODO better description of parameters
#' @param object a `cosme_model` object
#' @param ... additional arguments passed to `print`
#' @export
summary.cosme_fit <- function(object, ...) {
    object |> estimates() |> print(...)
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

    ## Get if info theory has been fit
    has_info <- !is.null(object$info)

    ## Get parameter estimates
    f <- parameterEstimates(object$freq)
    b <- parameterEstimates(object$bayes)

    ## Extract Bayesian posterior intervals
    # Have to capture output to prevent printing
    tmp <- capture.output(x <- summary(object$bayes), nullfile())
    # Extract intervals |> make numbers |> make matrix -> ci
    x[, c('pi.lower', 'pi.upper')] |> as.numeric() |> matrix(ncol = 2) -> bpi
    # Change NA to estimate
    is.na(bpi) |> ifelse(b$est, bpi) -> bpi
    # blavaan::blavInspect(object$bayes, 'hpd') # Doesnt match same rows columns

    ## Extract information theory
    if (has_info) {
        i <- vapply(object$info$mod_list[[1]], \(.) .@ParTable$est, numeric(NROW(f)))
        # intervals
        iint <- t(apply(i, 1, quantile, prob = c(0.025, 0.975, 0.5)))
    }

    ## Return
    structure(
        list(
            estimate = cbind(
                freq  = f$est,
                bayes = b$est,
                info  = if(has_info) iint[,3] else NULL
            ),
            interval_lower = cbind(
                freq  = f$ci.lower,
                bayes = bpi[,1],
                info  = if(has_info) iint[,1] else NULL
            ),
            interval_upper = cbind(
                freq  = f$ci.upper,
                bayes = bpi[,2],
                info  = if(has_info) iint[,2] else NULL
            )
        ),
        class = 'cosme_est',
        pname = f[,c('lhs', 'op', 'rhs')],
        ov    = c(
            lavNames(object$freq, type = 'ov.nox'),
            lavNames(object$freq, type = 'lv.nox')
        ),
        free  = is.na(object$freq@ParTable$ustart),
        has_info = has_info
    )
}

#' Print Estimates from `comse_est`
#' @description
#' A short description...
#'
# TODO better description of parameters
#' @param object a `cosme_model` object
#' @param ... Further arguments passed to or from other methods
#' @param nd  Number of decimal places to print (Defaults 3)
#' @import cli
#' @export
print.cosme_est <- function(x, ..., nd = 3L) {

    # Obtain parameter names and observed variables
    x |> attr('pname') -> pname
    x |> attr('ov') -> ov
    x |> attr('free') -> free

    # Obtain if it has info theory or not
    x |> attr('has_info') -> has_info

    # Obtain values
    (
        x
        |> do.call(what = cbind)
        |> round(nd)
        |> format(width = if(has_info) NULL else 10)
    ) -> values

    # Get max width
    values |> nchar() |> max() -> nwidth

    # Subset based on type
    latent <- which(pname$op == '=~')
    regression <- which(pname$op == '~')
    intercept  <- which(pname$op == '~1')
    covariance <- which(pname$op == '~~' & pname$lhs != pname$rhs)
    variance <- which(pname$op == '~~' & pname$lhs == pname$rhs)

    # Rearrange pname for duplicates
    for (i in seq_along(latent)) {
        if (i == length(latent)) break;
        for (j in (i+1):length(latent)) {
            si <- latent[i]
            sj <- latent[j]
            if (pname$lhs[si] == pname$lhs[sj]) pname$lhs[sj] <- ''
        }
    }
    for (i in seq_along(covariance)) {
        if (i == length(covariance)) break;
        for (j in (i+1):length(covariance)) {
            si <- covariance[i]
            sj <- covariance[j]
            if (pname$lhs[si] == pname$lhs[sj]) pname$lhs[sj] <- ''
        }
    }
    for (i in seq_along(regression)) {
        if (i == length(regression)) break;
        for (j in (i+1):length(regression)) {
            si <- regression[i]
            sj <- regression[j]
            if (pname$lhs[si] == pname$lhs[sj]) pname$lhs[sj] <- ''
        }
    }

    # Obtain padding
    padding <- max(max(nchar(pname$lhs) + 3L), 10L)
    ncol <- if(has_info) 3 else 2

    # Obtain printing functions
    funcs <- get_printfunc(nwidth, has_info)
    estimate_header1 <- funcs$estimate_header1
    estimate_header2 <- funcs$estimate_header2
    print_values     <- funcs$print_values

    # Begin printing
    cli::cli_h2('Estimates Summary')
    if (length(latent) > 0) {
        cli::cli_h3('Latent Variables:')
        cat(strrep(' ', padding + 4))
        cat(estimate_header1(), sep = '  ', fill = T)
        cat(strrep(' ', padding + 4))
        print_values(estimate_header2())
        for(i in latent) {
            val <- pname$lhs[i]
            val_r <- pname$rhs[i]
            if (nchar(val) != 0) {
                cat(formatC(paste0('  ', val, ' =~'), width = padding + 4, flag = '-'))
                cat(rep(strrep(' ', nwidth*ncol + ncol - 1), 3), sep = ' |', fill = T)
            }
            cat(formatC(paste0('    ', val_r), width = padding + 4, flag = '-'))
            print_values(values[i,], free = free[i])
        }
    }
    if (length(regression) > 0) {
        cli::cli_h3('Regressions:')
        cat(strrep(' ', padding + 4))
        cat(estimate_header1(), sep = '  ', fill = T)
        cat(strrep(' ', padding + 4))
        print_values(estimate_header2())
        for(i in regression) {
            val <- pname$lhs[i]
            val_r <- pname$rhs[i]
            if (nchar(val) != 0) {
                cat(formatC(paste0('  ', val, ' ~'), width = padding + 4, flag = '-'))
                cat(rep(strrep(' ', nwidth*ncol + ncol - 1), 3), sep = ' |', fill = T)
            }
            cat(formatC(paste0('    ', val_r), width = padding + 4, flag = '-'))
            print_values(values[i,], free = free[i])
        }
    }
    if (length(covariance) > 0) {
        cli::cli_h3('Covariances:')
        cat(strrep(' ', padding + 4))
        cat(estimate_header1(), sep = '  ', fill = T)
        cat(strrep(' ', padding + 4))
        print_values(estimate_header2())
        for(i in covariance) {
            val <- pname$lhs[i]
            val_r <- pname$rhs[i]
            if (nchar(val) != 0) {
                cat(formatC(paste0('  ', val, ' ~~'), width = padding + 4, flag = '-'))
                cat(rep(strrep(' ', nwidth*ncol + ncol - 1), 3), sep = ' |', fill = T)
            }
            cat(formatC(paste0('    ', val_r), width = padding + 4, flag = '-'))
            print_values(values[i,], free = free[i])
        }
    }
    if (length(intercept) > 0) {
        cli::cli_h3('Intercepts:')
        cat(strrep(' ', padding + 4))
        cat(estimate_header1(), sep = '  ', fill = T)
        cat(strrep(' ', padding + 4))
        print_values(estimate_header2())
        for(i in intercept) {
            if (pname$lhs[i] %in% ov) {
                cat(formatC(paste0('   .', pname$lhs[i]), width = padding + 4, flag = '-'))
            } else {
                cat(formatC(paste0('    ', pname$lhs[i]), width = padding + 4, flag = '-'))
            }
            print_values(values[i,], free = free[i])
        }
    }
    if (length(variance) > 0) {
        cli::cli_h3('Variances:')
        cat(strrep(' ', padding + 4))
        cat(estimate_header1(), sep = '  ', fill = T)
        cat(strrep(' ', padding + 4))
        print_values(estimate_header2())
        for(i in variance) {
            if (pname$lhs[i] %in% ov) {
                cat(formatC(paste0('   .', pname$lhs[i]), width = padding + 4, flag = '-'))
            } else {
                cat(formatC(paste0('    ', pname$lhs[i]), width = padding + 4, flag = '-'))
            }
            print_values(values[i,], free = free[i])
        }
    }

    invisible(x)
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
        object$bayes, baseline.model = attr(object, 'b_base'),
        fit.measures = c('BRMSEA', 'BCFI')
    )

    # Return output
    structure(
        list(
            # DOES NOT CURRENTLY WORK DUE TO BUG IN semTools
            # Need to wait for semTools to update
            # freq  = list(lavaan::fitmeasures(object$freq), semTools::moreFitIndices(object$freq)),
            freq  = lavaan::fitmeasures(
                object$freq, c(
                    "chisq", "df", "pvalue", "cfi", "rmsea",
                    "rmsea.ci.lower", "rmsea.ci.upper", "srmr"
                )
            ),
            bayes = c(
                ppp = lavaan::fitMeasures(object$bayes, c("ppp"))[[1]],
                cfi = mean(bfit@indices$BRMSEA),
                cfi = quantile(bfit@indices$BCFI, probs = 0.025),
                cfi = quantile(bfit@indices$BCFI, probs = 0.975),
                rmsea = mean(bfit@indices$BRMSEA),
                rmsea = quantile(bfit@indices$BRMSEA, probs = 0.025),
                rmsea = quantile(bfit@indices$BRMSEA, probs = 0.975)
            ),
            info  = if(is.null(object$info)) {
                NA

            } else {
                apply(object$info$fit_list[[1]], 2, quantile, probs = seq(0.1, 0.9, by = 0.1))
            }
        ),
        class = 'cosme_mfit'
    )
}
