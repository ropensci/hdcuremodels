#' Summarize a Fitted Mixture Cure Object.
#'
#' @description
#' \code{summary} method for a mixturecure object fit using \code{curegmifs},
#' \code{cureem}, \code{cv_curegmifs}, or \code{cv_cureem}.
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs},
#' \code{cureem}, \code{cv_curegmifs}, or \code{cv_cureem}.
#' @param ... other arguments.
#'
#' @return prints the number of non-zero coefficients from the incidence and
#' latency portions of the fitted mixture cure model when using the minimum AIC
#' to select the final model. When fitting a model using \code{curegmifs} or
#' \code{cureem} the summary function additionally prints results associated
#' with the following model selection methods: the step and value that maximizes
#' the log-likelihood; the step and value that minimizes the AIC, modified AIC
#' (mAIC), corrected AIC (cAIC), BIC, modified BIC (mBIC), and extended BIC
#' (EBIC).
#'
#' @export
#'
#' @seealso \code{\link{curegmifs}}, \code{\link{cureem}},
#' \code{\link{coef.mixturecure}}, \code{\link{plot.mixturecure}},
#' \code{\link{predict.mixturecure}}
#' @keywords methods
#' @method summary mixturecure
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
#' summary(fit)
summary.mixturecure <- function(object, ...) {
  if (!(class(object) %in% "mixturecure")) {
    stop("class of object must be mixturecure")
  }
  message("Mixture cure model fit using the ", object$method, " algorithm \n")
  if (object$cv == FALSE) {
    select <- select_model(object, model_select = "AIC")
    model_select_AIC <- which.min(select$AIC)
    model_select_cAIC <- which.min(select$cAIC)
    model_select_mAIC <- which.min(select$mAIC)
    model_select_BIC <- which.min(select$BIC)
    model_select_mBIC <- which.min(select$mBIC)
    model_select_EBIC <- which.min(select$EBIC)
    model_select_logLik <- which.max(select$logLik)
    message(
      "Number of non-zero incidence covariates at minimum AIC: ",
      sum(object$b_path[model_select_AIC,] != 0), "\n"
    )
    message(
      "Number of non-zero latency covariates at minimum AIC: ",
      sum(object$beta_path[model_select_AIC,] != 0), "\n"
    )
  if (object$cv == FALSE) {
    message("Optimal step for selected information criterion: ", object$method, " algorithm \n")
    message(
      "\t at step    = ", model_select_logLik, " logLik     = ",
      max(select$logLik), "\n"
    )
    message(
      "\t at step    = ", model_select_AIC, " AIC        = ",
      min(select$AIC), "\n"
    )
    message(
      "\t at step    = ", model_select_mAIC, " mAIC        = ",
      min(select$mAIC), "\n"
    )
    message(
      "\t at step    = ", model_select_cAIC, " cAIC        = ",
      min(select$cAIC), "\n"
    )
    message(
      "\t at step    = ", model_select_BIC, " BIC        = ",
      min(select$BIC), "\n"
    )
    message(
      "\t at step    = ", model_select_mBIC, " mBIC        = ",
      min(select$mBIC), "\n"
    )
    message(
      "\t at step    = ", model_select_EBIC, " EBIC        = ",
      min(select$EBIC), "\n"
    )
  }
  } else {
    message("using cross-validation \n")
    message(
      "Number of non-zero incidence covariates: ",
      sum(object$b != 0), "\n"
    )
    message(
      "Number of non-zero latency covariates: ",
      sum(object$beta != 0), "\n"
    )
  }
}
