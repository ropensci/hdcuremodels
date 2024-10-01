#' AUC for cure prediction using mean score imputation
#'
#' @description
#' This function calculates the AUC for cure prediction using the mean score
#' imputation (MSI) method proposed by Asano et al.
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs},
#' \code{cureem}, \code{cv_curegmifs}, \code{cv_cureem}.
#' @param newdata an optional data.frame that minimally includes the incidence
#' and/or latency variables to use for predicting the response. If omitted, the
#' training data are used.
#' @param cure_cutoff cutoff value for cure, used to produce a proxy for the
#' unobserved cure status; default is 5.
#' @param model_select a case-sensitive parameter for models fit using \code{curegmifs} or \code{cureem}
#' any step along the solution path can be selected. The default is
#' \code{model_select = "AIC"} which calculates the predicted values using
#' the coefficients from the model having the lowest AIC. Other options are
#' \code{model_select = "mAIC"} for the modified AIC,
#' \code{model_select = "cAIC"} for the corrected AIC,
#' \code{model_select = "BIC"}, \code{model_select = "mBIC"}
#' for the modified BIC, \code{model_select = "EBIC"} for the extended BIC,
#' \code{model_select = "logLik"} for the step that maximizes the
#' log-likelihood, or any numeric value from the solution path.
#' This option has no effect for objects fit using \code{cv_curegmifs} or
#' \code{cv_cureem}.
#'
#' @return Returns the AUC value for cure prediction using the mean score
#' imputation (MSI) method.
#' @export
#'
#' @references Asano, J., Hirakawa, H., Hamada, C. (2014) Assessing the
#' prediction accuracy of cure in the Cox proportional hazards cure model:
#' an application to breast cancer data. \emph{Pharmaceutical Statistics},
#' \bold{13}:357--363.
#'
#' @srrstats {G1.0} *Statistical Software should list at least one primary reference from published academic literature.*
#'
#' @seealso \code{\link{concordance_mcm}}
#'
#' @keywords univar
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' testing <- temp$testing
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'                   data = training, x_latency = training,
#'                   model = "weibull", thresh = 1e-4, maxit = 2000,
#'                   epsilon = 0.01, verbose = FALSE)
#' auc_mcm(fit, model_select = "cAIC")
#' auc_mcm(fit, newdata = testing)
auc_mcm <- function(object, newdata, cure_cutoff = 5, model_select = "AIC") {
  if (!("mixturecure" %in% class(object)))
    stop("class of object must be mixturecure")
  no_data <- (missing(newdata) || is.null(newdata))
  if (no_data) {
    testing_time <- object[["y"]][, 1]
    testing_delta <- object[["y"]][, 2]
    x_p <- object[["x_incidence"]]
  } else {
    df <- model.frame(parse(text = object$call, keep.source = FALSE)[[2]],
                      data = newdata)
    testing_time <- df[, 1][, 1]
    testing_delta <- df[, 1][, 2]
    x_p <- as.matrix(model.frame(as.formula(paste(" ~ ", paste(
              colnames(object$x_incidence), collapse = "+"))), newdata))
  }
  if (!is.null(x_p) && identical(x_p, object$x_incidence)) {
    if (object$scale) {
      sd <- apply(x_p, 2, sd)
      for (i in seq_len(dim(x_p)[2])) {
        if (sd[i] == 0) {
          x_p[, i] <- scale(x_p[, i], center = TRUE,
                            scale = FALSE)
        } else {
          x_p[, i] <- scale(x_p[, i], center = TRUE,
                            scale = TRUE)
        }
      }
    }
  } else if (!is.null(x_p) && object$scale) {
    newx <- rbind(object$x_incidence, x_p)
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
    x_p <- as.matrix(newx[-(1:dim(object$x_incidence)[1]), , drop = FALSE])
  }
  testing_n <- length(testing_time)
  v <- rep(0, testing_n)
  y <- rep(999, testing_n)
  y[testing_time > cure_cutoff] <- 0
  y[testing_time <= cure_cutoff & testing_delta == 1] <- 1
  v[y < 2] <- 1
  # Extract b_p_hat and intercept at model_select
  if (!object$cv) {
    if (is.character(model_select)) {
      model_select <- c("AIC", "BIC", "logLik", "cAIC", "mAIC",
                        "mBIC", "EBIC")[pmatch(model_select, c("AIC", "BIC",
                        "logLik", "cAIC", "mAIC", "mBIC", "EBIC"))]
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
      mAIC <- (2 + 2 * log(p/.5)) * df - 2 * logLik
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
    itct_hat <- object$b0_path[model_select]
    b_p_hat <- object$b_path[model_select, , drop = FALSE]

  } else {
    itct_hat <- object$b0
    b_p_hat <- matrix(object$b, ncol = 1)
  }
  if (all(b_p_hat == 0)) {
    xb <- rep(itct_hat, testing_n)
  }  else {
    xb <- itct_hat + x_p[, b_p_hat != 0, drop = FALSE] %*% b_p_hat[b_p_hat != 0]
  }
  p_hat <- 1 / (1 + exp(-xb))
  temp <- v * y + (1 - v) * p_hat
  temp1 <- temp[order(p_hat, decreasing = TRUE)]
  temp_f <- v * (1 - y) + (1 - v) * (1 - p_hat)
  temp1f <- temp_f[order(p_hat, decreasing = TRUE)]
  TPR <- c(0, cumsum(temp1) / cumsum(temp1)[testing_n])
  FPR <- c(0, cumsum(temp1f) / cumsum(temp1f)[testing_n])
  height <- (TPR[-1] + TPR[-length(TPR)]) / 2
  width <- diff(FPR)
  auc_msi <- sum(height * width)
  return(auc_msi)
}
