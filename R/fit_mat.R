# Internal function used for fitting models
# Modified from ockhamSEM (https://github.com/falkcarl/ockhamSEM)

#' @author Brian T. Keller, Carl F. Falk, and Michael Muthukrishna.
#' @description
#' This code has been modified from the original code in ockhamSEM
#' (\url{https://github.com/falkcarl/ockhamSEM})
#' @importFrom lavaan lavaan lavInspect fitMeasures
#' @noRd
fit_mat <- function(mat, lavmodel, vnames, saveModel, fit.measure) {
    # Add names to matrix
    colnames(mat) <- rownames(mat) <- vnames

    ## From original `fitmod` code (Carl F. Falk and Michael Muthukrishna)
    # added from ShortForm Tabu code; nuke all starting values? Does this work? Not with older lavaan version
    lavmodel@ParTable$est <- NULL
    lavmodel@ParTable$start <- NULL

    # replaced with update()?; might be slower, but more stable? (doesn't work with older lavaan version)
    # res<-try({my_fitted_model<-update(lavmodel, sample.cov=temp_matrix)

    ## Fit model with trycatch
    tryCatch(
        {
            lavmodel@Options$se <- "none"
            lavmodel@Options$start <- "default"
            lavmodel@Options$do.fit <- TRUE

            my_fitted_model <- lavaan(
                sample.cov = mat,
                sample.nobs = lavmodel@SampleStats@nobs |> unlist(),
                slotOptions = lavmodel@Options,
                slotParTable = lavmodel@ParTable,
                slotCache = lavmodel@Cache
            )

            fit_out <- if (my_fitted_model@optim$converged) {
                fitMeasures(my_fitted_model, fit.measure)
            } else {
                rep(NA, length(fit.measure))
            }

            # regardless, save fitted lavaan model for later inspection
            if (saveModel) {
                attr(fit_out, "model") <- my_fitted_model
            }
            # Return fit_out
            fit_out
        },
        error = function(e) rep(NA, length(fit.measure))
    )
}
