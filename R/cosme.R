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
    if (!is.model(model)) {
        throw_error(c(
            "The {.cls model} object is not properly constructed.",
            "x" = "The object is not of class {.cls cosme_model}"
        ))
    }
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
