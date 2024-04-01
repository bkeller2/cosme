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


#' Internal function for generating print functions
#' @noRd
get_printfunc <- function(nw, info) {
    force(nw); force(info)
    # Return functions
    list(
        estimate_header1 = function() {
            ncol <- if(info) 3 else 2
            c(' Estimate', ' Lower Bound (2.5%)', ' Upper Bound (97.5%)') |>
                format(width = nw*ncol + ncol - 1, justify = 'centre')
        },
        estimate_header2 = function() {
            rep(c('Freq', 'Bayes', if(info) 'Info' else NULL), 3) |>
                formatC(width = nw, flag = '+')
        },
        print_values = function(values, free = TRUE) {
            # Handle if info exists or not
            sel1 <- if(info) 1:3 else 1:2
            sel2 <- if(info) 4:6 else 3:4
            sel3 <- if(info) 7:9 else 5:6

            cat(values[sel1])
            cat(' |')
            if (free) cat(values[sel2])
            else cat(strrep(' ', nchar(values[sel2])))
            cat(' |')
            if (free) cat(values[sel3], fill = T)
            else cat(strrep(' ', nchar(values[sel3])), fill = T)
        }
    )
}

