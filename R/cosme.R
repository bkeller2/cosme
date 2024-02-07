#' Fit cosme model
#' @description
#' A short description...
#'
# TODO better description of parameters
#' @param model a `cosme_model` object
#' @param data a `data.frame`
#' @param ... passed onto other stuff
#' @import future
#' @export
cosme <- function(model, data, ...) {
    # Check that it is a model
    if (!is.model(model)) model <- make_model(model)

    # Run three fits
    # future evaluation of freq and bayes only
    freq  <- future(fit_freq(model, data, ...))
    bayes <- future(fit_bayes(model, data, ...))
    info  <- fit_info(model, data, ...)

    # Create `cosme_fit` and return
    structure(
        list2env(
            list(
                freq  = value(freq, stdout = FALSE),
                bayes = value(bayes, stdout = FALSE),
                info  = info
            ),
            parent = emptyenv()
        ),
        class = c("cosme_fit", "environment")
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
#' @import lavaan
#' @export
estimates <- function(object) {
    # Check that it is a `cosme_fit`
    if (!is.fit(object)) {
        throw_error(c(
            "x" = "The object is not of class {.cls cosme_fit}"
        ))
    }
    # Return output
    list(
        freq  = lavaan::parameterEstimates(object$freq),
        bayes = lavaan::parameterTable(object$bayes),
        info  = modlist_summary(object$info$mod_list)
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
    # Return output
    list(
        # DOES NOT CURRENTLY WORK DUE TO BUG IN semTools
        # Need to wait for semTools to update
        # freq  = list(lavaan::fitmeasures(object$freq), semTools::moreFitIndices(object$freq)),
        freq  = list(lavaan::fitmeasures(object$freq)),
        # Takes awhile
        bayes = list(lavaan::fitMeasures(object$bayes), blavaan::blavFitIndices(object$bayes)),
        info  = list(NULL) # TODO implement
    )
}
