#' Summarize a Fitted Mixture Cure Object.
#'
#' @description
#' \code{summary} method for a mixturecure object fit using \code{curegmifs}, \code{cureem}, \code{cv_curegmifs}, or \code{cv_cureem}.
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs}, \code{cureem}, \code{cv_curegmifs}, or \code{cv_cureem}.
#' @param ... other arguments.
#'
#' @return prints the following items extracted from the object fit using \code{curegmifs} or \code{cureem}: the step and value that maximizes the log-likelihood; the step and value that minimizes the AIC, modified AIC (mAIC), corrected AIC (cAIC), BIC, modified BIC (mBIC), and extended BIC (EBIC). Returns log-likelihood, AIC, and BIC if the object was fit using \code{cv_curegmifs} or \code{cv_cureem} at the optimal cross-validated values if no FDR control; the number of non-zero incidence and latency variables is returned when cross-validation is used together with FDR control.
#'
#' @export
#'
#' @seealso \code{\link{curegmifs}}, \code{\link{cureem}}, \code{\link{coef.mixturecure}}, \code{\link{plot.mixturecure}}, \code{\link{predict.mixturecure}}
#' @keywords methods
#' @method summary mixturecure
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(N = 100, J = 10, nTrue = 10, A = 1.8)
#' training <- temp$Training
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'                        data = training, x.latency = training,
#'                        model = "weibull", thresh = 1e-4, maxit = 2000,
#'                        epsilon = 0.01, verbose = FALSE)
#' summary(fit)
summary.mixturecure <-
  function (object, ...)
  {
    if (!(class(object) %in% "mixturecure")) stop("class of object must be mixturecure")
    if (!object$cv) {
       if (!is.null(object$x.incidence)) {
        vars.inc <- apply(object$b_path, 1, function(x) sum(x !=
                                                              0))
      }
      else {
        vars.inc <- 0
      }
      if (!is.null(object$x.latency)) {
        vars.lat <- apply(object$beta_path, 1, function(x) sum(x !=
                                                                 0))
      }
      else {
        vars.lat <- 0
      }
    }
    else {
      if (!is.null(object$x.latency)) {
        vars.lat <- sum(object$beta != 0)
      }
      else {
        vars.lat <- 0
      }
      if (!is.null(object$x.incidence)) {
        vars.inc <- sum(object$b != 0)
      }
      else {
        vars.inc <- 0
      }
    }
    if (object$model == "weibull") {
      df <- vars.inc + vars.lat + 3
    }
    else if (object$model == "exponential") {
      df <- vars.inc + vars.lat + 2
    }
    else if (object$model == "cox") {
      df <- vars.inc + vars.lat + 1
    }
    if (object$method == "EM")
      logLik <- object$logLik.inc + object$logLik.lat
    else logLik <- object$logLik
    p <- dim(object$x.incidence)[2] + dim(object$x.latency)[2]
    AIC <- 2 * df - 2 * logLik
    #cAIC <- AIC+(2*df^2+6*df+4)/(length(object$y)-df-2)
    cAIC<-AIC+(2*df*(df+1))/(length(object$y)-df-1) # https://www.mathworks.com/help/econ/information-criteria.html
    mAIC <- (2+2*log(p/.5)) * df - 2 * logLik
    BIC <- df * (log(length(object$y))) -  2 * logLik
    mBIC <- df * (log(length(object$y)) + 2*log(p/4)) - 2 * logLik
    EBIC <- log(length(object$y)) * df + 2*(1-.5)*log(choose(p, df)) - 2 * logLik
    if (object$cv == FALSE) {
      model.select.AIC <- which.min(AIC)
      model.select.cAIC <- which.min(cAIC)
      model.select.mAIC <- which.min(mAIC)
      model.select.BIC <- which.min(BIC)
      model.select.mBIC <- which.min(mBIC)
      model.select.EBIC <- which.min(EBIC)
      model.select.logLik = which.max(logLik)
    }
    cat("Mixture cure model fit using the", object$method, "algorithm \n")
    if (object$cv == FALSE) {
      cat("at step    = ", model.select.logLik, "logLik     = ", logLik[model.select.logLik], "\n")
      cat("at step    = ", model.select.AIC, "AIC        = ", AIC[model.select.AIC], "\n")
      cat("at step    = ", model.select.mAIC, "mAIC        = ", mAIC[model.select.mAIC], "\n")
      cat("at step    = ", model.select.cAIC, "cAIC        = ", cAIC[model.select.cAIC], "\n")
      cat("at step    = ", model.select.BIC, "BIC        = ", BIC[model.select.BIC], "\n")
      cat("at step    = ", model.select.mBIC, "mBIC        = ", mBIC[model.select.mBIC], "\n")
      cat("at step    = ", model.select.EBIC, "EBIC        = ", EBIC[model.select.EBIC], "\n")
      cat("\n")
    }
    else if (!object$fdr.control) {
      cat("using cross-validation \n")
      cat("logLik     = ", logLik, "\n")
      cat("AIC        = ", AIC, "\n")
      cat("BIC        = ", BIC, "\n")
      cat("\n")
    } else {
      cat("Number of non-zero incidence covariates", sum(object$b!=0), "\n")
      cat("Number of non-zero latency covariates", sum(object$beta!=0), "\n")
    }
  }
