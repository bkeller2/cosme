#' Generate comse model base class
#' @noRd
make_model <- function(model) {
    # TODO check if constructed appropriately

    # Return object
    structure(model, class = "cosme_model")
}

#' Check if `cosme_model`
#' @noRd
is.model <- function(x) {
    inherits(x, "cosme_model")
}

# Convert `cosme_model` to character
as.character.cosme_model <- function(x, ...) {
    as.character.default(x)
}

#' Convert `cosme_model` to `lavaan`
#' @noRd
as.lavaan <- function(x) {
    # Check that it is a model
    if (!is.model(x)) {
        throw_error(c(
            "The {.cls model} object is not properly constructed.",
            "x" = "The object is not of class {.cls cosme_model}"
        ))
    }

    # Return
    x |> as.character()
}
