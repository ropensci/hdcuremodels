#' Print the contents of a mixture cure fitted object
#'
#' @description
#' This function prints the first several incidence and latency coefficients,
#' the rate (when fitting an exponential or Weibull mixture cure model), and
#' alpha (when fitting a Weibull mixture cure model). This function returns
#' the fitted object invisible to the user.
#'
#' @param x a \code{mixturecure} object resulting from \code{curegmifs},
#' \code{cureem}, \code{cv_cureem}, or \code{cv_curegmifs}.
#' @param max maximum number of rows in a matrix or elements in a vector to
#' display
#' @param ... other arguments.
#'
#' @note The contents of a \code{mixturecure} fitted object differ depending
#' upon whether the EM (\code{cureem}) or GMIFS (\code{curegmifs}) algorithm is
#' used for model fitting or if cross-validation is used. Also, the output
#' differs depending upon whether \code{x_latency} is specified in the model
#' (i.e., variables are included in the latency portion of the model fit) or
#' only \code{terms} on the right hand side of the equation are included (i.e.,
#' variables are included in the incidence portion of the model).
#'
#' @seealso \code{\link{curegmifs}}, \code{\link{cureem}},
#' \code{\link{coef.mixturecure}}, \code{\link{summary.mixturecure}},
#' \code{\link{plot.mixturecure}}, \code{\link{predict.mixturecure}}
#'
#' @return prints coefficient estimates for the incidence portion of the model
#' and if included, prints the coefficient estimates for the latency portion of
#' the model. Also prints rate for exponential and Weibull models and scale
#' (alpha) for the Weibull mixture cure model. Returns all objects fit using
#' \code{cureem}, \code{curegmifs}, \code{cv_cureem}, or \code{cv_curegmifs}.
#' @export
#'
#' @importFrom utils head
#'
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {RE4.17} *Model objects returned by Regression Software should implement or appropriately extend a default `print` method which provides an on-screen summary of model (input) parameters and (output) coefficients.*
#' @keywords methods
#' @method print mixturecure
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
#' print(fit)
print.mixturecure <- function(x, max = 6, ...)  {
  cat("mixturecure object fit using ", x$model, x$method, "algorithm\n")
  cat("\n")
  if (x$cv) {
    cat("$b\n", head(x$b, n = max), "\n")
    cat(length(x$b) - max, "more elements\n\n")
    cat("$beta\n", head(x$beta, n = max), "\n")
    cat(length(x$beta) - max, "more elements\n\n")
    if (x$model %in% c("exponential", "weibull")) {
      cat("$rate\n", head(x$rate, n = max), "\n\n")
    }
    if (x$model %in% c("weibull")) {
      cat("$alpha\n", head(x$alpha, n = max), "\n")
    }
  } else {
    cat("$b_path\n")
    print.default(x$b_path[1 : max, 1 : max])
    cat(dim(x$b_path)[1] - max, "more rows\n")
    if (dim(x$b_path)[2] - max > 0)
      cat(dim(x$b_path)[2] - max, "more columns\n\n")
    cat("$beta_path\n")
    print.default(x$beta_path[1 : max, 1 : max])
    cat(dim(x$beta_path)[1] - max, "more rows\n")
    if (dim(x$beta_path)[2] - max > 0)
      cat(dim(x$beta_path)[2] - max, "more columns\n\n")
    if (x$model %in% c("exponential", "weibull")) {
      cat("$rate\n", head(x$rate, n = max), "\n")
      cat(length(x$rate) - max, "more elements\n\n")
    }
    if (x$model %in% c("weibull")) {
      cat("$alpha\n", head(x$alpha, n = max), "\n")
      cat(length(x$alpha) - max, "more elements\n\n")
    }
  }
  invisible(x)  # Ensure the function returns the object invisibly
}
