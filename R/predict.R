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
#' @export
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'                   data = training, x_latency = training,
#'                   model = "weibull", thresh = 1e-4, maxit = 2000,
#'                   epsilon = 0.01, verbose = FALSE)
#' predict_train <- predict(fit)
#' names(predict_train)
#' testing <- temp$testing
#' predict_test <- predict(fit, newdata = testing)
predict.mixturecure <- function(object, newdata, model_select = "AIC", ...) {
  if (!("mixturecure" %in% class(object)))
    stop("class of object must be mixturecure")
  no_data <- (missing(newdata) || is.null(newdata))
  if (no_data) {
    #y <- object$y
    x_inc <- object$x_incidence
    x_lat <- object$x_latency
  } else {
    #if (!is.null(newx_incidence))
    #    if (newx_incidence == ~1) {
    #        m <- model.frame(newx_incidence)
    #    }
    #    else {
    x_inc <- as.matrix(model.frame(as.formula(paste(" ~ ",
            paste(colnames(object$x_incidence), collapse = "+"))), newdata))
    #    }
    #if (!is.null(newx_latency))
    #    if (newx_latency == ~1) {
    #        m <- model.frame(newx_latency)
    #    }
    #    else {
    x_lat <- as.matrix(model.frame(as.formula(paste(" ~ ",
              paste(colnames(object$x_latency), collapse = "+"))), newdata))
    #    }
    #if (is.null(newx_incidence) & is.null(newx_latency)) {
    #    x_inc <- object$x_incidence
    #    x_lat <- object$x_latency
    #}
  }
  #n <- max(dim(x_inc)[1], dim(x_lat)[1])
  if (!is.null(x_inc) && identical(x_inc, object$x_incidence)) {
    if (object$scale) {
      sd <- apply(x_inc, 2, sd)
      for (i in seq_len(dim(x_inc)[2])) {
        if (sd[i] == 0) {
          x_inc[, i] <- scale(x_inc[, i], center = TRUE,
                              scale = FALSE)
        } else {
          x_inc[, i] <- scale(x_inc[, i], center = TRUE,
                              scale = TRUE)
        }
      }
    }
  } else if (!is.null(x_inc) && object$scale) {
    newx <- rbind(object$x_incidence, x_inc)
    sd <- apply(newx, 2, sd)
    for (i in seq_len(dim(newx)[2])) {
      if (sd[i] == 0) {
        newx[, i] <- scale(newx[, i], center = TRUE,
                           scale = FALSE)
      } else {
        newx[, i] <- scale(newx[, i], center = TRUE,
                           scale = TRUE)
      }
    }
    x_inc <- as.matrix(newx[-(seq_len(dim(object$x_incidence)[1])), ,
                            drop = FALSE])
  }
  if (!is.null(x_lat) && identical(x_lat, object$x_latency)) {
    if (object$scale) {
      sd <- apply(x_lat, 2, sd)
      for (i in seq_len(dim(x_lat)[2])) {
        if (sd[i] == 0) {
          x_lat[, i] <- scale(x_lat[, i], center = TRUE,
                              scale = FALSE)
        } else {
          x_lat[, i] <- scale(x_lat[, i], center = TRUE,
                              scale = TRUE)
        }
      }
    }
  } else if (!is.null(x_lat) && object$scale) {
    newx <- rbind(object$x_latency, x_lat)
    sd <- apply(newx, 2, sd)
    for (i in seq(dim(newx)[2])) {
      if (sd[i] == 0) {
        newx[, i] <- scale(newx[, i], center = TRUE,
                           scale = FALSE)
      } else {
        newx[, i] <- scale(newx[, i], center = TRUE,
                           scale = TRUE)
      }
    }
    x_lat <- as.matrix(newx[-(seq_len(dim(object$x_latency)[1])), ,
                            drop = FALSE])
  }
  if (!object$cv) {
    if (is.character(model_select)) {
      model_select <- c("AIC", "BIC", "logLik", "cAIC", "mAIC",
                        "mBIC", "EBIC")[pmatch(model_select,
                                               c("AIC", "BIC", "logLik", "cAIC",
                                                 "mAIC",
                                                 "mBIC", "EBIC"))]
      if (any(!model_select %in% c("AIC", "BIC", "logLik", "cAIC", "mAIC",
                                   "mBIC", "EBIC")))
        stop("model_select must be either 'AIC', 'BIC', 'logLik', 'cAIC',
             'mAIC', 'mBIC', or 'EBIC' ")
      if (object$method == "EM") {
        logLik <- object$logLik_inc + object$logLik_lat
      } else {
        logLik <- object$logLik
      }
      if (!is.null(object$x_incidence)) {
        vars_inc <- apply(object$b_path, 1, function(x) sum(x !=  0))
      } else {
        vars_inc <- 0
      }
      if (!is.null(object$x_latency)) {
        vars_lat <- apply(object$beta_path, 1, function(x) sum(x != 0))
      } else {
        vars_lat <- 0
      }
      if (object$model == "weibull") {
        df <- vars_inc + vars_lat + 3
      } else if (object$model == "exponential") {
        df <- vars_inc + vars_lat + 2
      } else if (object$model == "cox") {
        df <- vars_inc + vars_lat + 1
      }
      p <- dim(object$x_incidence)[2] + dim(object$x_latency)[2]
      AIC <- 2 * df - 2 * logLik
      cAIC <- AIC + (2 * df * (df + 1)) / (length(object$y) - df - 1)
      mAIC <- (2 + 2 * log(p / .5)) * df - 2 * logLik
      BIC <- df * (log(length(object$y))) -  2 * logLik
      mBIC <- df * (log(length(object$y)) + 2 * log(p / 4)) - 2 * logLik
      EBIC <- log(length(object$y)) * df + 2 * (1 - .5) * log(choose(p, df)) -
        2 * logLik
      if (model_select == "AIC") {
        model_select <- which.min(AIC)
      } else if (model_select == "BIC") {
        model_select <- which.min(BIC)
      } else if (model_select == "mAIC") {
        model_select <- which.min(mAIC)
      } else if (model_select == "mBIC") {
        model_select <- which.min(mBIC)
      } else if (model_select == "EBIC") {
        model_select <- which.min(EBIC)
      } else if (model_select == "cAIC") {
        model_select <- which.min(cAIC)
      } else if (model_select == "logLik") {
        model_select <- which.max(logLik)
      }
    }
    if (is.null(x_inc))
      p_hat <- 1 / (1 + exp(-object$b0_path[model_select]))
    else p_hat <- as.numeric(1 / (1 + exp(-object$b0_path[model_select] -
                                            x_inc %*%
                                              t(object$b_path[model_select, ,
                                                              drop = FALSE]))))
    if (is.null(x_lat)) {
      w_beta <- rep(0, dim(x_inc)[1])
    } else {
      w_beta <- as.numeric(x_lat %*% t(object$beta_path[model_select,
                                                        , drop = FALSE]))
      latency <- ifelse(w_beta < 0, "low risk", "high risk")
    }
    output <- list(p_uncured = p_hat, linear_latency = w_beta,
                   latency_risk = latency)
  } else {
    if (is.null(x_inc))
      p_hat <- 1 / (1 + exp(-object$b0))
    else p_hat <- 1 / (1 + exp(-object$b0 - x_inc %*% matrix(object$b,
                                                             ncol = 1)))
    if (is.null(x_lat)) {
      w_beta <- rep(0, dim(x_inc)[1])
    } else {
      w_beta <- as.numeric(x_lat %*% matrix(object$beta,
                                            ncol = 1))
      latency <- ifelse(w_beta < 0, "low risk", "high risk")
    }
    output <- list(p_uncured = as.numeric(p_hat), linear_latency = w_beta,
                   latency_risk = latency)
  }
  output
}
