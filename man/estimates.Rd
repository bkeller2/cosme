% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cosme.R
\name{estimates}
\alias{estimates}
\title{Obtain Estimates from \code{comse_fit}}
\usage{
estimates(object, std = FALSE)
}
\arguments{
\item{object}{a \code{cosme_model} object}

\item{std}{Logical. Standardized solution or not. See details}
}
\description{
A short description...
}
\details{
Returns lavaan's \code{std.all} method of standardization.
See \code{\link[lavaan:standardizedSolution]{lavaan::standardizedSolution}} for more information.
}
