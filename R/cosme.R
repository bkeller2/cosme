#' Fit cosme model
#' @description
#' A short description...
#'
# TODO better description of parameters
#' @param model a `cosme_model` object
#' @param data a `data.frame`
#' @param ... passed onto other stuff
#' @export
cosme <- function(model, data, ...) {
    # Check that it is a model
    if (!is.model(model)) {
        throw_error(c(
            "The {.cls model} object is not properly constructed.",
            "x" = "The object is not of class {.cls cosme_model}"
        ))
    }

    # Create `cosme_fit` and return
    structure(
        list2env(
            list(
                freq  = fit_freq(model, data, ...),
                bayes = fit_bayes(model, data, ...),
                info  = fit_info(model, data, ...)
            ),
            parent = emptyenv()
        ),
        class = "cosme_fit"
    )
}
