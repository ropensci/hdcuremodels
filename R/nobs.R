#' Number of observations for mixturecure objects
#'
#' @param object An object of class `mixturecure`.
#' @param ... other arguments.
#'
#' @return number of subjects in the dataset.
#' @export
#' @keywords methods
#' @method nobs mixturecure
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
#' nobs(fit)
nobs.mixturecure <- function(object, ...) {
  if (!("mixturecure" %in% class(object))) {
    stop("Error: class of object must be mixturecure")
  }
  if (is.null(object$y)) {
    d1 <- NULL
  }
  else d1 <- length(object$y)
  d1
}
