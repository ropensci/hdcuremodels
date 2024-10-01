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
#' @return prints the following items extracted from the object fit using
#' \code{curegmifs} or \code{cureem}: the step and value that maximizes the
#' log-likelihood; the step and value that minimizes the AIC, modified AIC
#' (mAIC), corrected AIC (cAIC), BIC, modified BIC (mBIC), and extended BIC
#' (EBIC). Returns log-likelihood, AIC, and BIC if the object was fit using
#' \code{cv_curegmifs} or \code{cv_cureem} at the optimal cross-validated values
#' if no FDR control; the number of non-zero incidence and latency variables is
#' returned when cross-validation is used together with FDR control.
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
#'                        data = training, x_latency = training,
#'                        model = "weibull", thresh = 1e-4, maxit = 2000,
#'                        epsilon = 0.01, verbose = FALSE)
#' summary(fit)
summary.mixturecure <- function(object, ...) {
  if (!(class(object) %in% "mixturecure"))
    stop("class of object must be mixturecure")
  if (!object$cv) {
    if (!is.null(object$x_incidence)) {
      vars_inc <- apply(object$b_path, 1, function(x) sum(x != 0))
    } else {
      vars_inc <- 0
    }
    if (!is.null(object$x_latency)) {
      vars_lat <- apply(object$beta_path, 1, function(x) sum(x != 0))
    } else {
      vars_lat <- 0
    }
  } else {
    if (!is.null(object$x_latency)) {
      vars_lat <- sum(object$beta != 0)
    } else {
      vars_lat <- 0
    }
    if (!is.null(object$x_incidence)) {
      vars_inc <- sum(object$b != 0)
    } else {
      vars_inc <- 0
    }
  }
  if (object$model == "weibull") {
    df <- vars_inc + vars_lat + 3
  } else if (object$model == "exponential") {
    df <- vars_inc + vars_lat + 2
  } else if (object$model == "cox") {
    df <- vars_inc + vars_lat + 1
  }
  if (object$method == "EM")
    logLik <- object$logLik_inc + object$logLik_lat
  else logLik <- object$logLik
  p <- dim(object$x_incidence)[2] + dim(object$x_latency)[2]
  AIC <- 2 * df - 2 * logLik
  cAIC <- AIC + (2 * df * (df + 1)) / (length(object$y) - df - 1)
  mAIC <- (2 + 2 * log(p / .5)) * df - 2 * logLik
  BIC <- df * (log(length(object$y))) -  2 * logLik
  mBIC <- df * (log(length(object$y)) + 2 * log(p / 4)) - 2 * logLik
  EBIC <- log(length(object$y)) * df + 2 * (1 - .5) * log(choose(p, df)) -
    2 * logLik
  if (object$cv == FALSE) {
    model_select_AIC <- which.min(AIC)
    model_select_cAIC <- which.min(cAIC)
    model_select_mAIC <- which.min(mAIC)
    model_select_BIC <- which.min(BIC)
    model_select_mBIC <- which.min(mBIC)
    model_select_EBIC <- which.min(EBIC)
    model_select_logLik <- which.max(logLik)
  }
  message("Mixture cure model fit using the ", object$method, " algorithm \n")
  if (object$cv == FALSE) {
    message("at step    = ", model_select_logLik, " logLik     = ",
            logLik[model_select_logLik], "\n")
    message("at step    = ", model_select_AIC, " AIC        = ",
            AIC[model_select_AIC], "\n")
    message("at step    = ", model_select_mAIC, " mAIC        = ",
            mAIC[model_select_mAIC], "\n")
    message("at step    = ", model_select_cAIC, " cAIC        = ",
            cAIC[model_select_cAIC], "\n")
    message("at step    = ", model_select_BIC, " BIC        = ",
            BIC[model_select_BIC], "\n")
    message("at step    = ", model_select_mBIC, " mBIC        = ",
            mBIC[model_select_mBIC], "\n")
    message("at step    = ", model_select_EBIC, " EBIC        = ",
            EBIC[model_select_EBIC], "\n")
    message("\n")
  } else if (!object$fdr_control) {
    message("using cross-validation \n")
    message("logLik     = ", logLik, "\n")
    message("AIC        = ", AIC, "\n")
    message("BIC        = ", BIC, "\n")
    message("\n")
  } else {
    message("Number of non-zero incidence covariates: ",
            sum(object$b != 0), "\n")
    message("Number of non-zero latency covariates: ",
            sum(object$beta != 0), "\n")
  }
}
