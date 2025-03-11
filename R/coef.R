#' Extract model coefficients from a fitted mixture cure object
#'
#' @description
#' \code{coef.mixturecure} is a generic function which extracts the model
#' coefficients from a fitted mixture cure model object fit using
#' \code{curegmifs}, \code{cureem}, \code{cv_curegmifs}, or \code{cv_cureem}.
#'
#' @aliases coefficients
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs},
#' \code{cureem}, \code{cv_curegmifs}, or \code{cv_cureem}.
#' @param model_select a case-sensitive parameter for models fit using \code{curegmifs} or \code{cureem}
#' any step along the solution path can be selected. The default is
#' \code{model_select = "AIC"} which calculates the predicted values using the
#' coefficients from the model having the lowest AIC. Other options are
#' \code{model_select = "mAIC"} for the modified AIC,
#' \code{model_select = "cAIC"} for the corrected AIC,
#' \code{model_select = "BIC"}, \code{model_select = "mBIC"} for the modified
#' BIC, \code{model_select = "EBIC"} for the extended BIC,
#' \code{model_select = "logLik"} for the step that maximizes the
#' log-likelihood, or any numeric value from the solution path. This option has
#' no effect for objects fit using \code{cv_curegmifs} or \code{cv_cureem}.
#' @param ... other arguments.
#'
#' @return a list of estimated parameters extracted from the model object using
#' the model selection criterion
#'
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {RE4.2} *Model coefficients (via `coef()` / `coefficients()`)*
#' @seealso \code{\link{curegmifs}}, \code{\link{cureem}},
#' \code{\link{summary.mixturecure}}, \code{\link{plot.mixturecure}},
#' \code{\link{predict.mixturecure}}
#' @export
#' @keywords methods
#' @method coef mixturecure
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'   data = training, x_latency = training,
#'   model = "weibull", thresh = 1e-4, maxit = 2000,
#'   epsilon = 0.01, verbose = FALSE
#' )
#' coef(fit)
coef.mixturecure <-
  function(object, model_select = "AIC", ...) {
    if (!(class(object) %in% "mixturecure")) {
      stop("Error: class of object must be mixturecure")
    }
    if (!object$cv) {
      if (!is.numeric(model_select))
        model_select <- select_model(object, model_select)$select
      if (is.null(object$x_latency) && !is.null(object$x_incidence)) {
        if (object$model == "exponential") {
          coef <- list(
            rate = object$rate_path[model_select],
            b0 = object$b0_path[model_select],
            beta_inc = object$b_path[model_select, ])
        }
        if (object$model == "weibull") {
          coef <- list(
            rate = object$rate_path[model_select],
            shape = object$alpha_path[model_select],
            b0 = object$b0_path[model_select],
            beta_inc = object$b_path[model_select, ])
        }
        if (object$model == "cox") {
          coef <- list(
            b0 = object$b0_path[model_select],
            beta_inc = object$b_path[model_select, ])
        }
      } else if (!is.null(object$x_latency) && is.null(object$x_incidence)) {
        if (object$model == "exponential") {
          coef <- list(
            rate = object$rate_path[model_select],
            beta_lat = object$beta_path[model_select, ]
          )
        }
        if (object$model == "weibull") {
          coef <- list(
            rate = object$rate_path[model_select],
            shape = object$alpha_path[model_select],
            beta_lat = object$beta_path[model_select, ]
          )
        }
        if (object$model == "cox") {
          coef <- list(beta_lat = object$beta_path[model_select, ])
        }
      } else if (!is.null(object$x_latency) && !is.null(object$x_incidence)) {
        if (object$model == "exponential") {
          coef <- list(
            rate = object$rate_path[model_select],
            b0 = object$b0_path[model_select],
            beta_inc = object$b_path[model_select, ],
            beta_lat = object$beta_path[model_select, ]
          )
        }
        if (object$model == "weibull") {
          coef <- list(
            rate = object$rate_path[model_select],
            shape = object$alpha_path[model_select],
            b0 = object$b0_path[model_select],
            beta_inc = object$b_path[model_select, ],
            beta_lat = object$beta_path[model_select, ]
          )
        }
        if (object$model == "cox") {
          coef <- list(
            b0 = object$b0_path[model_select],
            beta_inc = object$b_path[model_select, ],
            beta_lat = object$beta_path[model_select, ]
          )
        }
      }
      } else {
      if (object$model == "weibull") {
        if (!is.null(object$x_incidence) && !is.null(object$x_latency)) {
          coef <- list(
            rate = object$rate, shape = object$alpha,
            b0 = object$b0, beta_inc = object$b,
            beta_lat = object$beta
          )
        }
        if (is.null(object$x_incidence) && !is.null(object$x_latency)) {
          coef <- list(
            rate = object$rate, shape = object$alpha,
            beta_lat = object$beta
          )
        }
        if (!is.null(object$x_incidence) && is.null(object$x_latency)) {
          coef <- list(
            rate = object$rate, shape = object$alpha,
            b0 = object$b0, beta_inc = object$b
          )
        }
      } else if (object$model == "exponential") {
        if (!is.null(object$x_incidence) && !is.null(object$x_latency)) {
          coef <- list(
            rate = object$rate, b0 = object$b0,
            beta_inc = object$b, beta_lat = object$beta
          )
        }
        if (is.null(object$x_incidence) && !is.null(object$x_latency)) {
          coef <- list(rate = object$rate, beta_lat = object$beta)
        }
        if (!is.null(object$x_incidence) && is.null(object$x_latency)) {
          coef <- list(
            rate = object$rate, b0 = object$b0,
            beta_inc = object$b
          )
        }
      } else if (object$model == "cox") {
        if (!is.null(object$x_incidence) && !is.null(object$x_latency)) {
          coef <- list(
            b0 = object$b0, beta_inc = object$b,
            beta_lat = object$beta
          )
        }
        if (is.null(object$x_incidence) && !is.null(object$x_latency)) {
          coef <- list(beta_lat = object$beta)
        }
        if (!is.null(object$x_incidence) && is.null(object$x_latency)) {
          coef <- list(b0 = object$b0, beta_inc = object$b)
        }
      }
      }
    coef
  }
