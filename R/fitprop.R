# Main UI functions for ockhamSEM
# Taken from ockhamSEM (https://github.com/falkcarl/ockhamSEM)
# Modifications by Brian T. Keller

#' Run fit propensity analyses
#'
#' @param model Models of class lavaan for which the user would like to compare fit propensity.
#' @param fit.measure Character vector that indicates which fit measure to extract from fitted
#'   models. Possible options include anything returned by \code{\link[lavaan]{fitMeasures}} from the lavaan
#'   package when applied to fitted models.
#' @param rmethod String indicating the type of random correlation generation approach.
#'   Choices are \code{"mcmc"} (default), \code{"onion"}, and \code{"clustergen"}. See details.
#' @param reps Number of random correlation matrices to generate for fit propensity analysis.
#' @param onlypos Logical value indicating whether to generate correlation matrices. Note that if there
#'   are many variables, generation of correlation matrices and fitting models to them will be very, very
#'   computationally intensive.
#'   with positive manifold (\code{TRUE}); i.e., only positive relationships among variables.
#' @param mcmc.args Named list of arguments that controls options for
#'   \code{"mcmc"} correlation matrix generation. See details.
#' @param clustergen.args Named list of arguments that controls generation of
#'   correlation matrices if \code{"onion"} or \code{"clustergen"} is used. See details.
#' @param ... Passed to lavaan
#' @details Inspired by work by Preacher (2003, 2006) and Bonifay & Cai (2017),
#'   this function performs three steps for analyses to assess the fit propensity of competing
#'   structural equation models: 1. Randomly generate correlation (or covariance matrices);
#'   2. Fit models to each correlation matrix; and 3. Save a indices that could be used for
#'   evaluating model fit in subsequent summaries. Conceptually, models that exhibit better fit
#'   to such randomly generated data may have better fit propensity, and are therefore potentially
#'   less parsimonious.
#'
#'   Analyses are performed with the \code{lavaan} package, and fitted lavaan models of
#'   \code{\link[lavaan]{lavaan-class}} (e.g., created from
#'   \code{\link[lavaan]{cfa}}, \code{\link[lavaan]{sem}}, or \code{\link[lavaan]{lavaan}} functions)
#'   for the competing models must be passed as initial arguments
#'   to the function. Currently, only single-group models and those relying on ML estimation are
#'   supported. Otherwise, the underlying options for the fitted lavaan models will be re-used by
#'   the run.fitprop function for the fit propensity analyses. It is optional to save the randomly
#'   generated matrices from Step 1 and the models fit in Step 2. Follow-up summaries of results
#'   saved in Step 3 are provided by plot.fitprop and summary.fitprop functions.
#'
#'   Generation of random correlation matrices is provided using several approaches. The \code{"mcmc"}
#'   algorithm implements a Markov Chain Monte Carlo approach and was ported from Fortran code
#'   in Preacher (2003). For details on the algorithm's actual implementation, see Preacher (2003),
#'   Falk and Muthukrishna (in prep), or the source code for the mcmc function. If this algorithm
#'   is chosen, \code{mcmc.args} accepts a list that can modify some default settings. In particular,
#'   \code{iter} sets the total number of iterations to run (default = 5000000). If parallel processing
#'   is enabled, this number will be divided amonst the number of chains. \code{miniter} sets a
#'   minimum number of iterations per chain to avoid many processors leading to too few iterations per
#'   chain (default = 10000). \code{jmpsize} overrides the step size for each update to the candidate
#'   correlation matrix. Smaller step sizes typically lead to more acceptance and may be necessary for
#'   larger correlation matrices (default jump size depends on the number of variables). Though, in
#'   general the MCMC algorithm becomes more difficult to work well with many variables.
#'
#'   The \code{"onion"} method is one approach that relies on work of Joe (2006) and
#'   Lewandowski, Kurowick, and Joe (2009); matrices are generated recursively, one variable at a
#'   time. The onion method is computationally more efficient than the MCMC algorithm. Under the
#'   hood, the \code{\link[clusterGeneration]{genPositiveDefMat}} function in the clusterGeneration package is used, with default
#'   arguments of \code{covMethod="onion"}, \code{eta=1}, and \code{rangeVar=c(1,1)}. These arguments ensure that the
#'   Onion method is used, generation is uniform over the space of positive definite matrices
#'   (but see note on positive manifold below), and with unit variances.
#'
#'   An additional option \code{"clustergen"} is provided for direct interface with the \code{\link[clusterGeneration]{genPositiveDefMat}}
#'   function in the clusterGeneration package. A named list can be passed to \code{clustergen.args} to
#'   override any defaults used by \code{\link[clusterGeneration]{genPositiveDefMat}}, and the user is referred to documentation
#'   for that function. This allows, for example, generation using C-Vines, covariance matrices
#'   (i.e., variables that do not all have unit variances), and several other covaraince/correlation
#'   matrix generation techniques.
#'
#'   onlypos controls whether correlation matrices can have only positive correlations.
#'   The original MCMC algorith by Preacher (2003, 2006) generated correlation matrices with
#'   positive manifold only (i.e., only positive correlations). The algorithm is easily changed
#'   to allow also negative correlations. The Onion method and any functions from clusterGeneration
#'   by default generate matrices with both positive and negative correlations. To obtain
#'   matrices with positive manifold only, an ad-hoc correction is implemented for these latter
#'   approaches where the matrix is transformed: R = (R+1)/2. To our knowledge, there is no
#'   guarantee that this will result in uniform sampling from the space of all correlation matrices
#'   with positive manifold, yet fit propensity results for some examples are very similar to those
#'   of the MCMC algorithm.
#'
#' @slot fit_list A list of the same length as the number of models being compared. Each list
#'   contains a matrix with columns corresponding to each entry of \code{fit.measure} and for all replications.
#' @slot R A list of the same length as the number of replications, containing all correlation
#'   matrices that were used for the fit propensity analysis.
#' @slot mod_list A list of the same length as the number of models being compared. Each element
#'   contains a list of the same length as the number of replications and contains fitted lavaan
#'   models.
#' @references
#' Bonifay, W. E., & Cai, L. (2017). On the complexity of item response theory models. Multivariate Behavioral Research, 52(4), 465–484. \url{http://doi.org/10.1080/00273171.2017.1309262}
#'
#' Falk, C. F., & Muthukrishna, M. (in press). Parsimony in model selection: Tools for assessing fit propensity. Psychological Methods.
#'
#' Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random correlation matrices based on vines and extended onion method. Journal of Multivariate Analysis,100(9), 1989–2001. \url{http://doi.org/10.1016/j.jmva.2009.04.008}
#'
#' Joe, H. (2006). Generating random correlation matrices based on partial correlations. Journal of Multivariate Analysis, 97(10), 2177–2189. \url{http://doi.org/10.1016/j.jmva.2005.05.010}
#'
#' Preacher, K. J. (2003). The role of model complexity in the evaluation of structural equation models (PhD thesis). The Ohio State University.
#'
#' Preacher, K. J. (2006). Quantifying parsimony in structural equation modeling. Multivariate Behavioral Research, 41(3), 227–259. \url{http://doi.org/10.1207/s15327906mbr4103_1}
#'
#' Falk and Muthukrishna (2020) Parsimony in Model Selection: Tools for Assessing Fit Propensity.
#' @author Carl F. Falk and Michael Muthukrishna. MCMC ported from FORTRAN code by Kris Preacher (2003).
#'
#' Code modifications by Brian T. Keller
#' @examples
#' \donttest{
#' # Set up a covariance matrix to fit models to
#' p <- 3 # number of variables
#' temp_mat <- diag(p) # identity matrix
#' colnames(temp_mat) <- rownames(temp_mat) <- paste0("V", seq(1, p))
#'
#' # Define and fit two models using lavaan package
#' mod1a <- "V3 ~ V1 + V2
#'   V1 ~~ 0*V2"
#' mod2a <- "V3 ~ V1
#'   V2 ~ V3"
#'
#' mod1a.fit <- sem(mod1a, sample.cov = temp_mat, sample.nobs = 500)
#' mod2a.fit <- sem(mod2a, sample.cov = temp_mat, sample.nobs = 500)
#'
#' # Run fit propensity analysis a variety of different ways
#'
#' # Onion approach, only positive correlation matrices, save srmr
#' res <- run.fitprop(mod1a.fit, mod2a.fit,
#'     fit.measure = "srmr",
#'     rmethod = "onion", reps = 1000, onlypos = TRUE
#' )
#' summary(res)
#'
#' # Onion approach, save several fit indices
#' res <- run.fitprop(mod1a.fit, mod2a.fit,
#'     fit.measure = c("srmr", "cfi", "rmsea"),
#'     rmethod = "onion", reps = 1000
#' )
#' summary(res)
#' }
#' @return An object of class fitprop for which plot and summary methods are available. Some slots are listed in the Slots section.
#' @seealso \code{\link[ockhamSEM]{plot.fitprop}} \code{\link[ockhamSEM]{summary.fitprop}}
#' @importFrom lavaan lavInspect lavNames parTable lavaan
#' @import cli
#' @noRd
run_fitprop <- function(model, ...,
                        fit.measure = "srmr",
                        rmethod = c("onion", "mcmc", "clustergen"),
                        reps = 1000,
                        onlypos = FALSE,
                        mcmc.args = list(),
                        clustergen.args = list()) {
    fun.call <- match.call()

    if (class(model) == "lavaan") models <- list(model)
    else stop("Input model must be a fitted lavaan model")

    n.models <- length(models)

    for (lavmodel in models) {
        if (!class(lavmodel) == "lavaan") {
            stop("Input lavmodel must be a fitted lavaan model")
        }

        # Check for multiple group model
        if (lavInspect(lavmodel, "ngroups") > 1) {
            stop("Only single group models are currently supported")
        }
    }

    # Extract variable names
    vnames <- lavNames(models[[1]])

    # number of variables
    d <- length(vnames)

    # number of fit measures
    n.fit <- length(fit.measure)

    # check fit.measure
    if (!is.character(fit.measure)) {
        stop("fit.measure should be a character string or vector indicating the measure to be extracted from fitMeasures()")
    }

    # Set matrix for storing outputs
    fit_list <- vector("list", n.models)
    for (i in seq_along(fit_list)) {
        fit_list[[i]] <- matrix(0, reps, length(fit.measure))
        colnames(fit_list[[i]]) <- fit.measure
    }

    # Select appropriate control
    # Brian T. Keller
    control <- switch(
        rmethod,
        onion = onion.args.default(d, clustergen.args),
        mcmc = mcmc.args.default(1, d, mcmc.args),
        clustergen = clustergen.args.default(d, clustergen.args),
        stop("Unrecognized control option") # TODO change to cli
    )

    # Generate the matrices
    # Brian T. Keller
    out_mat <- lapply(
        seq_len(reps), genmat,
        rmethod = rmethod,
        control = control, onlypos = onlypos
    )

    # Do elements of @Data and @SampleStats depend on whether data is raw or cov matrix is used?
    # If so, pre-fit each model to cov matrix once
    # If not, this code isn't necessary
    # Update using lapply
    # Brian T. Keller
    models <- lapply(models, function(lavmodel) {
        temp_matrix <- lavInspect(lavmodel, "sampstat")$cov
        ptab <- parTable(lavmodel)
        opt <- lavmodel@Options
        opt$se <- "none"
        nobs <- lavInspect(lavmodel, "nobs")
        # update may not work with older lavaan versions
        # my_fitted_model<-update(lavmodel, sample.cov=temp_matrix, se="none")
        lavaan(
            model = ptab,
            sample.cov = temp_matrix,
            sample.nobs = nobs,
            slotOptions = opt
        )
    })

    # Fit model to each matrix
    # Update using lapply, future, and clean up code
    # Brian T. Keller
    cli::cli_alert_success("Matrices Generated.")
    cli::cli_alert('Fitting Models to Matrices.')

    mod_list <- vector("list", length(models))
    for (i in seq_along(models)) {
        lavmodel <- models[[i]]
        # with_progress({
        #     p <- progressor(along = out_mat)
        #     fit_tmp <- future.apply::future_lapply(
        #         out_mat,
        #         FUN = \(x, ...) {
        #             p(sprintf("x=%g", x))
        #             fit_mat(x, ...)
        #         },
        #         lavmodel = lavmodel,
        #         vnames = vnames,
        #         saveModel = TRUE,
        #         fit.measure = fit.measure,
        #         future.seed = FALSE
        #     )
        # }, handlers = handler_cli(format = paste0(
        #     "  Fitting {cli::pb_bar} {cli::pb_current}/{cli::pb_total} ",
        #     "| Elapsed: {cli::pb_elapsed}",
        #     "| ETA: {cli::pb_eta}"
        # ), show_after = F, clear = F))

        fit_tmp <- lapply(
            cli::cli_progress_along(
                seq_along(out_mat),
                clear = FALSE,
                format = paste0(
                    "  Fitting {cli::pb_bar} {cli::pb_current}/{cli::pb_total} ",
                    "| Elapsed: {cli::pb_elapsed}",
                    "| ETA: {cli::pb_eta}"
                )
            ),
            FUN = \(x, ...) {
                flush.console()
                fit_mat(out_mat[[x]], ...)
            },
            lavmodel = lavmodel,
            vnames = vnames,
            saveModel = TRUE,
            fit.measure = fit.measure
        )
        fit_list[[i]] <- matrix(
            unlist(fit_tmp),
            ncol = length(fit.measure),
            byrow = TRUE
        )
        colnames(fit_list[[i]]) <- fit.measure
        mod_list[[i]] <- lapply(fit_tmp, attr, "mod")
    }
    cli::cli_alert_success("Models fit to all matrices.")
    flush.console()

    # Process number of missing values (if any) per model and fit index
    na_list <- lapply(fit_list, is.na)
    complete_mat <- do.call(cbind, lapply(na_list, ifelse, yes = NA, no = 1))

    # Return class
    structure(list(
        fun.call = fun.call,
        reps = reps,
        rmethod = rmethod,
        onlypos = onlypos,
        nfit = n.fit,
        origmodels = models,
        na_list = na_list,
        complete_mat = complete_mat,
        fit_list = fit_list,
        R = out_mat,
        mod_list = mod_list
    ), class = "fitprop")
}

# #' Plot function for fitprop objects
# #' @param x Object of class fitprop, as created by run.fitprop function.
# #' @param ... Does nothing but to hopefully make this generic function pass R CMD check
# #' @param type What type of plot to produce? Options include \code{"ecdf"} or \code{"euler"}. \code{"nVennR"}
# #'   is not yet operational.
# #' @param whichmod Index number corresponding to which model(s) to include on the plot. Defaults to all models.
# #' @param whichfit Character vector indicating which indices of model fit to include. Defaults to all saved indices.
# #' @param savePlot Logical value indicating whether to save plot to a list (TRUE) or just produce plot to output.
# #' @param xlim Numeric vector of length 2 indicating the limits for the x-axis of the plot.
# #' @param samereps Logical value indicating whether to use only results from replications in which all selected models yielded results.
# #' @param cutoff Numeric vector indicating what cut value of the fit indice(s) to use for euler plots.
# #' @param lower.tail Logical vector indicating whether lower values of each fit index corresponds to good fit.
# #' @param mod.lab Optional character vector of labels for each model.
# #' @param mod.brewer.pal Optional string corresponding to the palette from RColorBrewer to use for the different models. e.g.,
# #'   see \code{\link[RColorBrewer]{display.brewer.all}}.
# #' @examples
# #' \donttest{
# #'
# #' library(ggplot2)
# #' # Set up a covariance matrix to fit models to
# #' p <- 3 # number of variables
# #' temp_mat <- diag(p) # identity matrix
# #' colnames(temp_mat) <- rownames(temp_mat) <- paste0("V", seq(1, p))
# #'
# #' # Define and fit two models using lavaan package
# #' mod1a <- "V3 ~ V1 + V2
# #'   V1 ~~ 0*V2"
# #' mod2a <- "V3 ~ V1
# #'   V2 ~ V3"
# #'
# #' mod1a.fit <- sem(mod1a, sample.cov = temp_mat, sample.nobs = 500)
# #' mod2a.fit <- sem(mod2a, sample.cov = temp_mat, sample.nobs = 500)
# #'
# #' # Run fit propensity analysis
# #' # Onion approach, save srmr and CFI
# #' res <- run.fitprop(mod1a.fit, mod2a.fit,
# #'     fit.measure = c("srmr", "cfi"),
# #'     rmethod = "onion", reps = 1000
# #' )
# #' summary(res)
# #'
# #' # Generate a variety of plots
# #' # ecdf - only srmr
# #' plot(res, whichfit = "srmr")
# #'
# #' # ecdf - both srmr and cfi, properly indicating that higher values of cfi are better
# #' plot(res, lower.tail = c(TRUE, FALSE))
# #'
# #' # ecdf change color palette, save plot to object, then display or change
# #' myplot <- plot(res,
# #'     whichfit = "cfi",
# #'     lower.tail = FALSE,
# #'     mod.brewer.pal = "Dark2",
# #'     savePlot = TRUE
# #' )
# #' myplot[[1]] # saved to the first slot
# #' myplot[[1]] + theme_bw() # change something, like the theme
# #'
# #' # euler plot
# #' # cfi - which models have cfi of .5 or better?
# #' plot(res, type = "euler", whichfit = "cfi", lower.tail = FALSE, cutoff = .5)
# #' }
# #' @author Carl F. Falk and Michael Muthukrishna. MCMC ported from FORTRAN code by Kris Preacher (2003).
# #' @references Falk and Muthukrishna (2020) Parsimony in Model Selection: Tools for Assessing Fit Propensity.
# #' @noRd
# #' @import ggplot2
# #' @importFrom graphics plot
# #' @importFrom stats na.omit
# #' @importFrom tidyr gather
# #' @importFrom eulerr euler
# #' @importFrom rlang .data
# plot.fitprop <- function(x, ..., type = c("ecdf", "euler", "nVennR"), whichmod = NULL, whichfit = colnames(x$fit_list[[1]]), savePlot = FALSE,
#                          xlim = c(0, 1), samereps = TRUE, cutoff = rep(.1, length(whichfit)), lower.tail = rep(TRUE, length(whichfit)),
#                          mod.lab = NULL, mod.brewer.pal = "Set1") {
#     type <- match.arg(type)
#
#     data <- x$fit_list
#     nmod <- length(data) # number of models
#     nfit <- ncol(data[[1]]) # number of fit measures
#     nrep <- nrow(data[[1]]) # number of available replications
#
#     if (is.null(whichfit)) {
#         whichfit <- colnames(data[[1]])
#     }
#     if (is.null(whichmod)) {
#         whichmod <- 1:nmod
#     }
#     if (is.null(mod.lab)) {
#         mod.lab <- paste0("Model ", whichmod)
#     }
#
#     plots <- list()
#     # loop over fit indices
#     m <- 1
#     for (fm in whichfit) {
#         # extract and format data
#         dat <- matrix(nrow = nrep, ncol = nmod)
#         j <- 1
#         for (mod in 1:nmod) {
#             dat[, j] <- data[[mod]][, fm]
#             j <- j + 1
#         }
#         dat <- as.data.frame(dat)
#         colnames(dat) <- mod.lab
#         dat$id <- 1:nrow(dat)
#
#         # generate plots
#         if (type == "ecdf") {
#             if (samereps) {
#                 dat <- na.omit(dat)
#             }
#             dat <- gather(dat, "variable", "value", whichmod)
#             graph <- ggplot(dat, aes(x = .data$value)) +
#                 stat_ecdf(aes(linetype = .data$variable, color = .data$variable), na.rm = TRUE, size = .7) +
#                 scale_color_brewer(palette = mod.brewer.pal)
#
#             graph <- graph + guides(color = guide_legend(title = "Model"))
#             graph <- graph + guides(linetype = guide_legend(title = "Model"))
#
#             if (lower.tail[m]) {
#                 graph <- graph + xlim(xlim[1], xlim[2])
#             } else {
#                 graph <- graph + xlim(xlim[2], xlim[1])
#             }
#             graph <- graph + xlab(fm) + ylab("Cumulative Probability") # theme(legend.title=element_blank())+
#         } else if (type == "euler") {
#             if (lower.tail[m]) {
#                 dat[, whichmod] <- dat[, whichmod] < cutoff # how many meet cutoff criterion
#             } else {
#                 dat[, whichmod] <- dat[, whichmod] > cutoff # how many meet cutoff criterion
#             }
#             tmp <- cbind(data.frame(Total = rep(TRUE, nrow(dat))), dat[, whichmod])
#             if (samereps) {
#                 tmp <- na.omit(tmp)
#             }
#             colnames(tmp) <- c("Total", whichmod)
#             eulerfit <- euler(tmp)
#             graph <- plot(eulerfit)
#         } else if (type == "nVennR") {
#             stop("nVennR is not yet functional")
#
#             dat <- na.omit(dat)
#             tmp <- list()
#             indx <- 1
#             for (j in whichmod) {
#                 if (lower.tail[m]) {
#                     tmp[[indx]] <- dat$id[dat[, j] < cutoff]
#                 } else {
#                     tmp[[indx]] <- dat$id[dat[, j] > cutoff]
#                 }
#                 indx <- indx + 1
#             }
#             tmp[[indx]] <- dat$id # last group is the total
#             # tmp<-cbind(data.frame(Total=rep(TRUE,nrow(dat))),dat[,whichmod])
#             # if(samereps){
#             #  tmp<-na.omit(tmp)
#             # }
#             # colnames(tmp)<-c("Total",whichmod)
#             # eulerfit<-euler(tmp)
#             # graph<-plot(eulerfit)
#             # graph<-plotVenn(tmp,nCycles=5000,showPlot=F,sNames=c(mod.lab,"Total"))
#         }
#
#         # do something with plot
#         if (savePlot) {
#             plots[[m]] <- graph
#         } else {
#             if (m > 1) {
#                 invisible(readline(prompt = "Press [enter] to continue"))
#             }
#             print(graph)
#         }
#         m <- m + 1
#     }
#     if (savePlot) {
#         invisible(plots)
#     }
# }

#' @author Carl F. Falk and Michael Muthukrishna. MCMC ported from FORTRAN code by Kris Preacher (2003).
#' @references Falk and Muthukrishna (2020) Parsimony in Model Selection: Tools for Assessing Fit Propensity.
#' @noRd
print.fitprop <- function(x, ...) {
    data <- x$fit_list
    nmod <- length(data) # number of models
    nfit <- ncol(data[[1]]) # number of fit measures
    nrep <- x$reps # number of available replications

    cat(
        "\n",
        # "Function call =           ", deparse(x$fun.call), "\n",
        "Number of fitted models = ", nmod, "\n",
        "Number of fit measures =  ", nfit, "\n",
        "Fit measures =            ", colnames(data[[1]]), "\n",
        "Number of replications =   ", nrep, "\n",
        "Options for R generation\n",
        "  method =                ", x$rmethod, "\n",
        "  Only positive R?        ", x$onlypos, "\n"
    )
}
