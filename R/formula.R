#' Extract model formula for mixturecure object
#'
#' @description Extract the model formula for `mixturecure` object
#' @param x an object from class `mixturecure`.
#' @param ... other arguments.
#'
#' @returns a formula representing the incidence and variables for the latency
#' portion of the model
#' @keywords methods
#' @method formula mixturecure
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @export
#'
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
#' formula(fit)
formula.mixturecure <- function(x, ...) {
  if (!("mixturecure" %in% class(x))) {
    stop("Error: class of object must be mixturecure")
  }
  if (is.null(x$call))
    form <- NULL
  else
    form <- formula(x$call)
  latency <- as.list(x$call)$x_latency
  incidence <- list(formula = form, data = as.list(x$call)$data)
  latency <- list(formula = form, data = as.list(x$call)$x_latency)
  list(incidence = incidence, latency = latency)
}
