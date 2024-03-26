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

#' Internal function to produce baseline model
#' @noRd
baseline_model <- function(model) {
    model |>
        lav_partable_independence() |>
        with({
            c(paste(lhs, op, rhs), paste(lhs, '~ 1'))
        }) |>
        make_model()
}

#' Internal function to obtain estimates header1
#' @noRd
estimate_header1 <- function(nw) {
    c(' Estimate', ' Lower Bound (2.5%)', ' Upper Bound (97.5%)') |>
        format(width = nw*3 + 2, justify = 'centre')
}
#' Internal function to obtain estimates header2
#' @noRd
estimate_header2 <- function(nw) {
    rep(c('Freq', 'Bayes', 'Info'), 3)  |> formatC(width = nw, flag = '+')
}

#' Internal function to print out values
#' @noRd
print_values <- function(values, free = TRUE) {
    cat(values[1:3])
    cat(' |')
    if (free) cat(values[4:6])
    else cat(strrep(' ', nchar(values[4:6])))
    cat(' |')
    if (free) cat(values[7:9], fill = T)
    else cat(strrep(' ', nchar(values[7:9])), fill = T)
}
