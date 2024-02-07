#' Internal error function
#' Wrapper for `cli_abort` to not specify the call
#' @noRd
throw_error <- function(message, ..., .envir = parent.frame(), .frame = .envir) {
    cli::cli_abort(message, ..., .envir = .envir, .frame = .frame, call = NULL)
}


#' Internal function for producing average estimates
#' @noRd
modlist_summary <- function(object) {
    # Make sure only one model compared
    if (length(object) != 1) {
        throw_error(c(
            "!" = "Incorrect number of models specified."
        ))
    }
    object <- object[[1]]

    # Make sure that there is more than one replication
    if (length(object) == 0) {
        throw_error(c(
            "x" = "Not enough replications were specified."
        ))
    }
    val <- sapply(object, \(.) parameterEstimates(.)$est)
    out <- cbind(
        mean = rowMeans(val),
        sd = apply(val, 1, sd),
        t(apply(val, 1, quantile))
    )

    # Set collapsed parameter names as row names
    rownames(out) <- apply(
        parameterEstimates(object[[1]])[,c('lhs', 'op', 'rhs')],
        1, paste0, collapse = ' '
    )
    return(out)
}
