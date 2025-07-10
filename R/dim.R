#' Dimension method for mixturecure objects
#'
#' @description Dimension method for \code{mixturecure} objects.
#' @param x An object of class `mixturecure`.
#'
#' @return \item{nobs}{ number of subjects in the dataset.}
#' @return \item{p_incidence}{ number of variables in the incidence portion of the model.}
#' @return \item{p_latency}{ number of variables in the latency portion of the model.}
#' @export
#' @keywords methods
#' @method dim mixturecure
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
#' dim(fit)
dim.mixturecure <- function(x) {
  if (!("mixturecure" %in% class(x))) {
    stop("Error: class of object must be mixturecure")
  }
  d1name <- "nobs"
  d2name <- "p_incidence"
  d3name <- "p_latency"
  if (is.null(x$y)) {
    d1 <- d1name <- NULL
  }
  else d1 <- length(x$y)
  if (is.null(x$x_incidence)) {
    d2 <- d2name <- NULL
  } else {
    d2 <- ncol(x$x_incidence)
  }
  if (is.null(x$x_latency)) {
    d3 <- d3name <- NULL
  } else {
    d3 <- ncol(x$x_latency)
  }
  dd <- c(d1, d2, d3)
  names(dd) <- c(d1name, d2name, d3name)
  dd
}
