#' Predicted probabilities for susceptibles, linear predictor for latency, and risk class for latency for mixture cure fit
#'
#' @description
#' This function returns a list the includes the predicted probabilities for susceptibles as well as the linear predictor for the latency distribution and a dichotomous risk for latency for a \code{curegmifs}, \code{cureem}, \code{cv_curegmifs} or \code{cv_cureem} fitted object.
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs}, \code{cureem}, \code{cv_curegmifs}, \code{cv_cureem}.
#' @param newdata an optional data.frame that minimally includes the incidence and/or latency variables to use for predicting the response. If omitted, the training data are used.
#' @param model.select for models fit using \code{curegmifs} or \code{cureem} any step along the solution path can be selected. The default is \code{model.select = "AIC"} which calculates the predicted values using the coefficients from the model having the lowest AIC. Other options are \code{model.select = "mAIC"} for the modified AIC,  \code{model.select = "cAIC"} for the corrected AIC, \code{model.select = "BIC"}, \code{model.select = "mBIC"} for the modified BIC, \code{model.select = "EBIC"} for the extended BIC, \code{model.select = "logLik"} for the step that maximizes the log-likelihood, or any numeric value from the solution path. This option has no effect for objects fit using \code{cv_curegmifs} or \code{cv_cureem}.
#' @param ... other arguments
#'
#' @return \item{p.uncured}{ a vector of probabilities from the incidence portion of the fitted model representing the P(uncured).}
#' @return \item{linear.latency}{ a vector for the linear predictor from the latency portion of the model.}
#' @return \item{latency.risk}{ a dichotomous class representing low (below the median) versus high risk for the latency portion of the model.}
#'
#' @seealso \code{\link{curegmifs}}, \code{\link{cureem}}, \code{\link{coef.mixturecure}}, \code{\link{summary.mixturecure}}, \code{\link{plot.mixturecure}}
#' @keywords methods
#' @method predict mixturecure
#'
#' @export
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(N = 100, J = 10, nTrue = 10, A = 1.8)
#' training <- temp$Training
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'                   data = training, x.latency = training,
#'                   model = "weibull", thresh = 1e-4, maxit = 2000,
#'                   epsilon = 0.01, verbose = FALSE)
#' predict.train <- predict(fit)
#' names(predict.train)
#' testing <- temp$Testing
#' predict.test <- predict(fit, newdata = testing)
predict.mixturecure <- function (object, newdata, model.select = "AIC", ...) {
  if (!("mixturecure" %in% class(object))) stop("class of object must be mixturecure")
  noData <- (missing(newdata) || is.null(newdata))
  if (noData) {
    y <- object$y
    x.inc <- object$x.incidence
    x.lat <- object$x.latency
  } else {
    #if (!is.null(newx.incidence))
    #    if (newx.incidence == ~1) {
    #        m <- model.frame(newx.incidence)
    #    }
    #    else {
    x.inc <- as.matrix(model.frame(as.formula(paste(" ~ ", paste(colnames(object$x.incidence), collapse= "+"))), newdata))
    #    }
    #if (!is.null(newx.latency))
    #    if (newx.latency == ~1) {
    #        m <- model.frame(newx.latency)
    #    }
    #    else {
    x.lat <- as.matrix(model.frame(as.formula(paste(" ~ ", paste(colnames(object$x.latency), collapse= "+"))), newdata))
    #    }
    #if (is.null(newx.incidence) & is.null(newx.latency)) {
    #    x.inc <- object$x.incidence
    #    x.lat <- object$x.latency
    #}
  }
  n <- max(dim(x.inc)[1], dim(x.lat)[1])
  if (!is.null(x.inc) & identical(x.inc, object$x.incidence)) {
    if (object$scale) {
      sd <- apply(x.inc, 2, sd)
      for (i in 1:dim(x.inc)[2]) {
        if (sd[i] == 0) {
          x.inc[, i] <- scale(x.inc[, i], center = TRUE,
                              scale = FALSE)
        }
        else {
          x.inc[, i] <- scale(x.inc[, i], center = TRUE,
                              scale = TRUE)
        }
      }
    }
  }
  else if (!is.null(x.inc) && object$scale) {
    newx <- rbind(object$x.incidence, x.inc)
    sd <- apply(newx, 2, sd)
    for (i in 1:dim(newx)[2]) {
      if (sd[i] == 0) {
        newx[, i] <- scale(newx[, i], center = TRUE,
                           scale = FALSE)
      }
      else {
        newx[, i] <- scale(newx[, i], center = TRUE,
                           scale = TRUE)
      }
    }
    x.inc <- as.matrix(newx[-(1:dim(object$x.incidence)[1]),,drop=FALSE])
  }
  if (!is.null(x.lat) & identical(x.lat, object$x.latency)) {
    if (object$scale) {
      sd <- apply(x.lat, 2, sd)
      for (i in 1:dim(x.lat)[2]) {
        if (sd[i] == 0) {
          x.lat[, i] <- scale(x.lat[, i], center = TRUE,
                              scale = FALSE)
        }
        else {
          x.lat[, i] <- scale(x.lat[, i], center = TRUE,
                              scale = TRUE)
        }
      }
    }
  }
  else if (!is.null(x.lat) && object$scale) {
    newx <- rbind(object$x.latency, x.lat)
    sd <- apply(newx, 2, sd)
    for (i in 1:dim(newx)[2]) {
      if (sd[i] == 0) {
        newx[, i] <- scale(newx[, i], center = TRUE,
                           scale = FALSE)
      }
      else {
        newx[, i] <- scale(newx[, i], center = TRUE,
                           scale = TRUE)
      }
    }
    x.lat <- as.matrix(newx[-(1:dim(object$x.latency)[1]),,drop=FALSE])
  }
  if (!object$cv) {
    if (is.character(model.select)) {
      model.select <- c("AIC", "BIC", "logLik", "cAIC", "mAIC", "mBIC", "EBIC")[pmatch(model.select,
                                                               c("AIC", "BIC", "logLik", "cAIC", "mAIC", "mBIC", "EBIC"))]
      if (any(!model.select%in%c("AIC", "BIC", "logLik", "cAIC", "mAIC", "mBIC", "EBIC")))
        stop("model.select must be either 'AIC', 'BIC', 'logLik', 'cAIC', 'mAIC', 'mBIC', or 'EBIC' ")
      if (object$method == "EM") {
        logLik <- object$logLik.inc + object$logLik.lat
      } else {
        logLik <- object$logLik
      }
      if (!is.null(object$x.incidence)) {
        vars.inc <- apply(object$b_path, 1, function(x) sum(x !=  0))
      } else {
        vars.inc <- 0
      }
      if (!is.null(object$x.latency)) {
        vars.lat <- apply(object$beta_path, 1, function(x) sum(x != 0))
      } else {
        vars.lat <- 0
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
      p <- dim(object$x.incidence)[2] + dim(object$x.latency)[2]
      AIC <- 2 * df - 2 * logLik
      #cAIC <- AIC+(2*df^2+6*df+4)/(length(object$y)-df-2)
      cAIC<-AIC+(2*df*(df+1))/(length(object$y)-df-1) # https://www.mathworks.com/help/econ/information-criteria.html
      mAIC <- (2+2*log(p/.5)) * df - 2 * logLik
      BIC <- df * (log(length(object$y))) -  2 * logLik
      mBIC <- df * (log(length(object$y)) + 2*log(p/4)) - 2 * logLik
      EBIC <- log(length(object$y)) * df + 2*(1-.5)*log(choose(p, df)) - 2 * logLik
      if (model.select == "AIC") {
        model.select = which.min(AIC)
      }
      else if (model.select == "BIC") {
        model.select = which.min(BIC)
      }
      else if (model.select == "mAIC") {
        model.select = which.min(mAIC)
      }
      else if (model.select == "mBIC") {
        model.select = which.min(mBIC)
      }
      else if (model.select == "EBIC") {
        model.select = which.min(EBIC)
      }
      else if (model.select == "cAIC") {
        model.select = which.min(cAIC)
      }
      else if (model.select == "logLik") {
        model.select = which.max(logLik)
      }
    }
    if (is.null(x.inc))
      p_hat = 1/(1 + exp(-object$b0_path[model.select]))
    else p_hat = as.numeric(1/(1 + exp(-object$b0_path[model.select] -
                                         x.inc %*% t(object$b_path[model.select, , drop = FALSE]))))
    if (is.null(x.lat)) {
      W_beta = rep(0, dim(x.inc)[1])
    }
    else {
      W_beta <- as.numeric(x.lat %*% t(object$beta_path[model.select,
                                                        , drop = FALSE]))
      latency <- ifelse(W_beta < 0, "low risk", "high risk")
    }
    output <- list(p.uncured = p_hat, linear.latency = W_beta,
                   latency.risk = latency)
  }
  else {
    if (is.null(x.inc))
      p_hat = 1/(1 + exp(-object$b0))
    else p_hat = 1/(1 + exp(-object$b0 - x.inc %*% matrix(object$b,
                                                          ncol = 1)))
    if (is.null(x.lat)) {
      W_beta = rep(0, dim(x.inc)[1])
    }
    else {
      W_beta <- as.numeric(x.lat %*% matrix(object$beta,
                                            ncol = 1))
      latency <- ifelse(W_beta < 0, "low risk", "high risk")
    }
    output <- list(p.uncured = as.numeric(p_hat), linear.latency = W_beta,
                   latency.risk = latency)
  }
  output
}
