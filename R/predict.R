#' Predicted probabilities for susceptibles, linear predictor for latency, and
#' risk class for latency for mixture cure fit
#'
#' @description
#' This function returns a list the includes the predicted probabilities for
#' susceptibles as well as the linear predictor for the latency distribution
#' and a dichotomous risk for latency for a \code{curegmifs}, \code{cureem},
#' \code{cv_curegmifs} or \code{cv_cureem} fitted object.
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs},
#' \code{cureem}, \code{cv_curegmifs}, \code{cv_cureem}.
#' @param newdata an optional data.frame that minimally includes the incidence
#' and/or latency variables to use for predicting the response. If omitted, the
#' training data are used.
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
#' @param ... other arguments
#'
#' @return \item{p_uncured}{ a vector of probabilities from the incidence
#' portion of the fitted model representing the P(uncured).}
#' @return \item{linear_latency}{ a vector for the linear predictor from the
#' latency portion of the model.}
#' @return \item{latency_risk}{ a dichotomous class representing low (below the
#' median) versus high risk for the latency portion of the model.}
#'
#' @seealso \code{\link{curegmifs}}, \code{\link{cureem}},
#' \code{\link{coef.mixturecure}}, \code{\link{summary.mixturecure}},
#' \code{\link{plot.mixturecure}}
#' @keywords methods
#' @method predict mixturecure
#'
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {RE4.9} *Modelled values of response variables.*
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
#' predict_train <- predict(fit)
#' names(predict_train)
#' testing <- temp$testing
#' predict_test <- predict(fit, newdata = testing)
predict.mixturecure <- function(object, newdata, model_select = "AIC", ...) {
  if (!("mixturecure" %in% class(object))) {
    stop("Error: class of object must be mixturecure")
  }
  no_data <- (missing(newdata) || is.null(newdata))
  if (no_data) {
    x_inc <- object$x_incidence
    x_lat <- object$x_latency
  } else {
    x_inc <- as.matrix(model.frame(as.formula(paste(
      " ~ ",
      paste(colnames(object$x_incidence), collapse = "+")
    )), newdata))
    x_lat <- as.matrix(model.frame(as.formula(paste(
      " ~ ",
      paste(colnames(object$x_latency), collapse = "+")
    )), newdata))
  }
  # Scale x_inc if scale is TRUE; match to train if test data used
  if (!is.null(x_inc) && identical(x_inc, object$x_incidence)) {
    if (object$scale) {
      x_inc <- self_scale(x_inc, object$scale)
    }
  } else if (!is.null(x_inc) && object$scale) {
    newx <- rbind(object$x_incidence, x_inc)
    newx <- self_scale(newx, object$scale)
    x_inc <- as.matrix(newx[-(seq_len(dim(object$x_incidence)[1])), ,
                          drop = FALSE
    ])
  }
  # Scale x_lat if scale is TRUE; match to train if test data used
  if (!is.null(x_lat) && identical(x_lat, object$x_latency)) {
    if (object$scale) {
      x_lat <- self_scale(x_lat, object$scale)
    }
  } else if (!is.null(x_lat) && object$scale) {
    newx <- rbind(object$x_latency, x_lat)
    newx <- self_scale(newx, object$scale)
    x_lat <- as.matrix(newx[-(seq_len(dim(object$x_latency)[1])), , drop = FALSE])
  }
  if (!object$cv) {
    if (!is.numeric(model_select))
      model_select <- select_model(object, model_select)$select
    if (is.null(x_inc)) {
      p_hat <- 1 / (1 + exp(-object$b0_path[model_select]))
    } else {
      p_hat <- as.numeric(1 / (1 + exp(-object$b0_path[model_select] -
        x_inc %*%
        t(object$b_path[model_select, ,
          drop = FALSE]))))
    }
    if (is.null(x_lat)) {
      w_beta <- rep(0, dim(x_inc)[1])
    } else {
      w_beta <- as.numeric(x_lat %*% t(object$beta_path[model_select, ,
        drop = FALSE]))
      latency <- ifelse(w_beta < 0, "low risk", "high risk")
    }
    output <- list(p_uncured = p_hat, linear_latency = w_beta,
      latency_risk = latency)
    } else {
      if (is.null(x_inc)) {
        p_hat <- 1 / (1 + exp(-object$b0))
      } else {
        p_hat <- 1 / (1 + exp(-object$b0 - x_inc %*% matrix(object$b,
        ncol = 1)))
      }
      if (is.null(x_lat)) {
        w_beta <- rep(0, dim(x_inc)[1])
      } else {
        w_beta <- as.numeric(x_lat %*% matrix(object$beta,
        ncol = 1))
        latency <- ifelse(w_beta < 0, "low risk", "high risk")
      }
      output <- list(
      p_uncured = as.numeric(p_hat), linear_latency = w_beta,
      latency_risk = latency)
    }
  output
}
