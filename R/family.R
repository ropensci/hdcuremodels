#' Return model family and fitting algorithm for mixturecure model fits
#'
#' @description Return model family and fitting algorithm for\code{mixturecure}
#' model fits.
#' @param object an object of class \code{mixturecure}
#' @param  ... other arguments.
#'
#' @returns the parametric or semi-parametric model fit and the fitting
#' algorithm.
#' @export
#' @method family mixturecure
#' @keywords methods
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @examples
#' library(survival)
#' withr::local_seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'   data = training, x_latency = training,
#'   model = "weibull", thresh = 1e-4, maxit = 2000,
#'   epsilon = 0.01, verbose = FALSE
#' )
#' family(fit)
family.mixturecure <- function(object, ...) {
  if (!("mixturecure" %in% class(object))) {
    stop("Error: class of object must be mixturecure")
  }
  model <- paste0(toupper(substr(object$model, 1, 1)), substr(object$model, 2, nchar(object$model)))
  cat("\n")
  cat("Family:", model, "\n")
  cat("Algorithm:", object$method, "\n\n")
}
