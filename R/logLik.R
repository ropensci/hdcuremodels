#' Log-likelihood for fitted mixture cure model
#'
#' @description
#' This function returns the log-likelihood for a user-specified model criterion
#' or step for a \code{curegmifs}, \code{cureem},
#' \code{cv_curegmifs} or \code{cv_cureem} fitted object.
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs},
#' \code{cureem}, \code{cv_curegmifs}, \code{cv_cureem}.
#' @param model_select either a case-sensitive parameter for models fit using
#' \code{curegmifs} or \code{cureem} or any numeric step along the solution path
#' can be selected. The default is \code{model_select = "AIC"} which calculates
#' the predicted values using the coefficients from the model having the lowest
#' AIC. The complete list of options are:
#' \itemize{
#'     \item \code{"AIC"} for the minimum AIC (default).
#'     \item \code{"mAIC"} for the minimum modified AIC.
#'     \item \code{"cAIC"} for the minimum corrected AIC.
#'     \item \code{"BIC"}, for the minimum BIC.
#'     \item \code{"mBIC"} for the minimum modified BIC.
#'     \item \code{"EBIC"} for the minimum extended BIC.
#'     \item \code{"logLik"} for the step that maximizes the
#' log-likelihood.
#'     \item \code{n} where n is any numeric value from the
#'     solution path.
#'   }
#' This option has no effect for objects fit using \code{cv_curegmifs} or
#' \code{cv_cureem}.
#' @param ... other arguments.
#'
#' @return log-likelihood of the fitted mixture cure model using the specified
#' criteria.
#'
#' @keywords methods
#' @method logLik mixturecure
#'
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
#' logLik(fit, model_select = "AIC")
logLik.mixturecure <- function(object, model_select = "AIC", ...) {
  if (!("mixturecure" %in% class(object))) {
    stop("Error: class of object must be mixturecure")
  }
  if (object$cv) {
    if (object$method == "GMIFS") {
      logLik <- object$logLik
    } else {
    logLik <- object$logLik.inc + object$logLik.lat
    }
  } else {
    if (is.numeric(model_select) && model_select > length(object$logLik))
      stop("Error: model_select step must be less than or equal to ", length(object$logLik))
    if (!is.numeric(model_select))
      model_select <- select_model(object, model_select)$select
    logLik <- object$logLik[model_select]
  }
  logLik
}
