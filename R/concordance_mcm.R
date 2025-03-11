#' C-statistic for mixture cure models
#'
#' @description
#' This function calculates the C-statistic using the cure status weighting
#' (CSW) method proposed by Asano and Hirakawa.
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs},
#'  \code{cureem}, \code{cv_curegmifs}, \code{cv_cureem}.
#' @param newdata an optional data.frame that minimally includes the incidence
#' and/or latency variables to use for predicting the response. If omitted, the
#' training data are used.
#' @param cure_cutoff cutoff value for cure, used to produce a proxy for the
#' unobserved cure status; default is 5.
#' @param model_select a case-sensitive parameter for models fit using
#' \code{curegmifs} or \code{cureem}
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
#'
#' @return value of C-statistic for the cure models.
#' @export
#'
#' @references Asano, J. and Hirakawa, H. (2017) Assessing the prediction
#' accuracy of a cure model for censored survival data with long-term survivors:
#' Application to breast cancer data. \emph{Journal of Biopharmaceutical
#' Statistics}, \bold{27}:6, 918--932.
#'
#' @srrstats {G1.0} *Statistical Software should list at least one primary reference from published academic literature.*
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {G1.3} *All statistical terminology should be clarified and unambiguously defined.*
#' @srrstats {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @seealso \code{\link{auc_mcm}}
#'
#' @keywords univar
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' testing <- temp$testing
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'   data = training, x_latency = training,
#'   model = "weibull", thresh = 1e-4, maxit = 2000,
#'   epsilon = 0.01, verbose = FALSE
#' )
#' concordance_mcm(fit, model_select = "cAIC")
#' concordance_mcm(fit, newdata = testing, model_select = "cAIC")
concordance_mcm <- function(object, newdata, cure_cutoff = 5,
  model_select = "AIC") {
  if (!("mixturecure" %in% class(object))) {
    stop("Error: class of object must be mixturecure")
  }
  no_data <- (missing(newdata) || is.null(newdata))
  if (no_data) {
    testing_time <- object[["y"]][, 1]
    testing_delta <- object[["y"]][, 2]
    x_p <- object[["x_incidence"]]
    w_p <- object[["x_latency"]]
  } else {
    df <- model.frame(parse(text = object$call, keep.source = FALSE)[[2]],
      data = newdata
    )
    testing_time <- df[, 1][, 1]
    testing_delta <- df[, 1][, 2]
    x_p <- as.matrix(model.frame(as.formula(paste(
      " ~ ",
      paste(colnames(object$x_incidence), collapse = "+")
    )), newdata))
    w_p <- as.matrix(model.frame(as.formula(paste(
      " ~ ",
      paste(colnames(object$x_latency), collapse = "+")
    )), newdata))
  }
  # Scale x_p if scale is TRUE; match to train if test data used
  if (!is.null(x_p) && identical(x_p, object$x_incidence)) {
    if (object$scale) {
      x_p <- self_scale(x_p, object$scale)
    }
  } else if (!is.null(x_p) && object$scale) {
    newx <- rbind(object$x_incidence, x_p)
    newx <- self_scale(newx, object$scale)
    x_p <- as.matrix(newx[-(seq_len(dim(object$x_incidence)[1])), ,
                          drop = FALSE
    ])
  }
  # Scale w_p if scale is TRUE; match to train if test data used
  if (!is.null(w_p) && identical(w_p, object$x_latency)) {
    if (object$scale) {
      w_p <- self_scale(w_p, object$scale)
    }
  } else if (!is.null(w_p) && object$scale) {
    newx <- rbind(object$x_latency, w_p)
    newx <- self_scale(newx, object$scale)
    w_p <- as.matrix(newx[-(seq_len(dim(object$x_latency)[1])), , drop = FALSE])
  }
  c_csw_num <- 0
  c_csw_denom <- 0
  testing_n <- length(testing_time)
  v <- rep(0, testing_n)
  y <- rep(999, testing_n)
  y[testing_time > cure_cutoff] <- 0
  y[testing_time <= cure_cutoff & testing_delta == 1] <- 1
  v[y < 2] <- 1
  # Extract b_p_hat and intercept at model_select when CV not used
  if (!object$cv) {
    if (!is.numeric(model_select))
      model_select <- select_model(object, model_select)$select
    itct_hat <- object$b0_path[model_select]
    b_p_hat <- object$b_path[model_select, , drop = FALSE]
    beta_p_hat <- object$beta_path[model_select, , drop = FALSE]
  } else {
    itct_hat <- object$b0
    b_p_hat <- matrix(object$b, ncol = 1)
    beta_p_hat <- matrix(object$beta, ncol = 1)
  }
  if (all(b_p_hat == 0)) {
    p_hat <- 1 / (1 + exp(-itct_hat))
  } else {
    p_hat <- 1 / (1 + exp(-itct_hat - x_p[, b_p_hat != 0, drop = FALSE] %*%
      b_p_hat[b_p_hat != 0]))
  }
  temp <- v * y + (1 - v) * as.vector(p_hat)
  if (all(beta_p_hat == 0)) {
    w_beta <- rep(0, testing_n)
  } else {
    w_beta <- w_p[, beta_p_hat != 0, drop = FALSE] %*%
      beta_p_hat[beta_p_hat != 0]
  }
  for (i in 1:testing_n) {
    for (j in 1:testing_n) {
      if (j == i || !testing_delta[i] || testing_time[i] > testing_time[j]) next
      I_ij <- testing_time[i] < testing_time[j] |
        (testing_time[i] == testing_time[j] & !testing_delta[j])
      if (!I_ij) next
      if (w_beta[i] > w_beta[j]) c_csw_num <- c_csw_num + temp[j]
      c_csw_denom <- c_csw_denom + temp[j]
    }
  }
  return(c_csw_num / c_csw_denom)
}
