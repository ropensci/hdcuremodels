#' C-statistic for mixture cure models
#'
#' @description
#' This function calculates the C-statistic using the cure status weighting (CSW) method proposed by Asano and Hirakawa.
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs}, \code{cureem}, \code{cv_curegmifs}, \code{cv_cureem}.
#' @param newdata an optional data.frame that minimally includes the incidence and/or latency variables to use for predicting the response. If omitted, the training data are used.
#' @param cure_cutoff cutoff value for cure, used to produce a proxy for the unobserved cure status; default is 5.
#' @param model.select for models fit using \code{curegmifs} or \code{cureem} any step along the solution path can be selected. The default is \code{model.select = "AIC"} which calculates the predicted values using the coefficients from the model having the lowest AIC. Other options are \code{model.select = "mAIC"} for the modified AIC,  \code{model.select = "cAIC"} for the corrected AIC, \code{model.select = "BIC"}, \code{model.select = "mBIC"} for the modified BIC, \code{model.select = "EBIC"} for the extended BIC, \code{model.select = "logLik"} for the step that maximizes the log-likelihood, or any numeric value from the solution path. This option has no effect for objects fit using \code{cv_curegmifs} or \code{cv_cureem}.
#'
#' @return value of C-statistic for the cure models.
#' @export
#'
#' @references Asano, J. and Hirakawa, H. (2017) Assessing the prediction accuracy of a cure model for censored survival data with long-term survivors: Application to breast cancer data. \emph{Journal of Biopharmaceutical Statistics}, \bold{27}:6, 918--932.
#'
#' @seealso \code{\link{AUC}}
#'
#' @keywords univar
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(N = 100, J = 10, nTrue = 10, A = 1.8)
#' training <- temp$Training
#' testing <- temp$Testing
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'                   data = training, x.latency = training,
#'                   model = "weibull", thresh = 1e-4, maxit = 2000,
#'                   epsilon = 0.01, verbose = FALSE)
#' concordance_mcm(fit)
#' concordance_mcm(fit, newdata = testing)
concordance_mcm <- function(object, newdata, cure_cutoff = 5, model.select = "AIC"){
  if (!("mixturecure" %in% class(object))) stop("class of object must be mixturecure")
  noData <- (missing(newdata) || is.null(newdata))
  if (noData) {
    testing_time <- object[["y"]][,1]
    testing_delta <- object[["y"]][,2]
    X_p <- object[["x.incidence"]]
    W_p <- object[["x.latency"]]
  } else {
    df <- model.frame(parse(text=object$call, keep.source=FALSE)[[2]], data=newdata)
    testing_time <- df[,1][,1]
    testing_delta <- df[,1][,2]
    X_p <- as.matrix(model.frame(as.formula(paste(" ~ ", paste(colnames(object$x.incidence), collapse= "+"))), newdata))
    W_p <- as.matrix(model.frame(as.formula(paste(" ~ ", paste(colnames(object$x.latency), collapse= "+"))), newdata))
  }
  if (!is.null(X_p) & identical(X_p, object$x.incidence)) {
    if (object$scale) {
      sd <- apply(X_p, 2, sd)
      for (i in 1:dim(X_p)[2]) {
        if (sd[i] == 0) {
          X_p[, i] <- scale(X_p[, i], center = TRUE,
                            scale = FALSE)
        }
        else {
          X_p[, i] <- scale(X_p[, i], center = TRUE,
                            scale = TRUE)
        }
      }
    }
  }
  else if (!is.null(X_p) && object$scale) {
    newx <- rbind(object$x.incidence, X_p)
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
    X_p <- as.matrix(newx[-(1:dim(object$x.incidence)[1]),,drop=FALSE])
  }
  if (!is.null(W_p) & identical(W_p, object$x.latency)) {
    if (object$scale) {
      sd <- apply(W_p, 2, sd)
      for (i in 1:dim(W_p)[2]) {
        if (sd[i] == 0) {
          W_p[, i] <- scale(W_p[, i], center = TRUE,
                            scale = FALSE)
        }
        else {
          W_p[, i] <- scale(W_p[, i], center = TRUE,
                            scale = TRUE)
        }
      }
    }
  }
  else if (!is.null(W_p) && object$scale) {
    newx <- rbind(object$x.latency, W_p)
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
    W_p <- as.matrix(newx[-(1:dim(object$x.latency)[1]),,drop=FALSE])
  }
  C_csw_num = 0
  C_csw_denom = 0
  testing_n = length(testing_time)
  v = rep(0, testing_n)
  y = rep(999, testing_n)
  y[testing_time>cure_cutoff] = 0
  y[testing_time<=cure_cutoff & testing_delta==1] = 1
  v[y<2] = 1
  # Extract b_p_hat and intercept at model.select
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
    itct_hat = object$b0_path[model.select]
    b_p_hat = object$b_path[model.select, , drop = FALSE]
    beta_p_hat = object$beta_path[model.select, , drop = FALSE]
  } else {
    itct_hat = object$b0
    b_p_hat = matrix(object$b, ncol = 1)
    beta_p_hat = matrix(object$beta, ncol = 1)
  }
  if(all(b_p_hat==0)) {
    p_hat = 1/(1+exp(-itct_hat))
  } else{
    p_hat = 1/(1+exp(-itct_hat- X_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0]))
  }
  temp = v*y + (1-v)*as.vector(p_hat)
  if(all(beta_p_hat==0)) {
     W_beta = rep(0, testing_n)
  } else{
    W_beta = W_p[,beta_p_hat!=0,drop=FALSE] %*% beta_p_hat[beta_p_hat!=0]
  }
  for(i in 1:testing_n)
    for(j in 1:testing_n){
      if (j==i | !testing_delta[i] | testing_time[i]>testing_time[j]) next
      I_ij = testing_time[i]<testing_time[j] | (testing_time[i]==testing_time[j] & !testing_delta[j])
      if (!I_ij) next
      if (W_beta[i]>W_beta[j]) C_csw_num = C_csw_num + temp[j]
      C_csw_denom = C_csw_denom + temp[j]
    }
  return(C_csw_num / C_csw_denom)
}
