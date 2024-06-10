#' Fit a penalized parametric mixture cure model using the GMIFS algorithm with cross-validation for model selection
#'
#' @description
#' Fits a penalized Weibull or exponential mixture cure model using the generalized monotone incremental forward stagewise (GMIFS) algorithm with k-fold cross-validation to select the optimal iteration step along the solution path. When FDR controlled variable selection is used, the model-X knockoffs method is applied and indices of selected variables are returned.
#'
#' @param formula an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The response must be a survival object as returned by the \code{Surv} function while the variables on the right side of the formula are the covariates that are included in the incidence portion of the model.
#' @param data a data.frame in which to interpret the variables named in the \code{formula} or in the \code{subset} argument.
#' @param subset an optional expression indicating which subset of observations to be used in the fitting process, either a numeric or factor variable should be used in subset, not a character variable. All observations are included by default.
#' @param x.latency specifies the variables to be included in the latency portion of the model and can be either a matrix of predictors, a model formula with the right hand side specifying the latency variables, or the same data.frame passed to the \code{data} parameter. Note that when using the model formula syntax for \code{x.latency} it cannot handle \code{x.latency = ~ .}.
#' @param model type of regression model to use for the latency portion of mixture cure model. Can be "weibull" or "exponential"; default is "weibull".
#' @param penalty.factor.inc vector of binary indicators representing the penalty to apply to each incidence coefficient: 0 implies no shrinkage and 1 implies shrinkage. If not supplied, 1 is applied to all incidence variables.
#' @param penalty.factor.lat vector of binary indicators representing the penalty to apply to each latency coefficient: 0 implies no shrinkage and 1 implies shrinkage. If not supplied, 1 is applied to all latency variables.
#' @param fdr.control logical, if TRUE, model-X knockoffs are used for FDR-controlled variable selection and indices of selected variables are returned (default is FALSE).
#' @param fdr numeric value in (0, 1) range specifying the target FDR level to use for variable selection when \code{fdr.control=TRUE} (default is 0.2).
#' @param epsilon small numeric value reflecting incremental value used to update a coefficient at a given step (default is 0.001).
#' @param thresh small numeric value. The iterative process stops when the differences between successive expected penalized complete-data log-likelihoods for both incidence and latency components are less than this specified level of tolerance (default is 10^-5).
#' @param scale logical, if TRUE the predictors are centered and scaled.
#' @param maxit integer specifying the maximum number of steps to run in the iterative algorithm (default is 10^4).
#' @param inits an optional list specifiying the initial value for the incidence intercept (\code{itct}), a numeric vector for the unpenalized incidence coefficients (\code{b_u}), and a numeric vector for unpenalized latency coefficients (\code{beta_u}), a numeric value for the rate parameter (\code{lambda}), and a numeric value for the shape parameter (\code{alpha}) when \code{model = "weibull"}. If not supplied or improperly supplied, initialization is automatically provided by the function.
#' @param n_folds an integer specifying the number of folds for the k-fold cross-valiation procedure (default is 5).
#' @param measure.inc character string specifying the evaluation criterion used in selecting the optimal \eqn{\lambda_b}. Can be "c" or "auc"; default is "c". If \code{measure.inc="c"}, the C-statistic using the cure status weighting (CSW) method proposed by Asano and Hirakawa (2017) is used to select both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}. If \code{measure.inc="auc"}, the AUC for cure prediction using the mean score imputation (MSI) method proposed by Asano et al. (2014) is used to select \eqn{\lambda_b} while the C-statistic with CSW is used for \eqn{\lambda_{\beta}}.
#' @param one.se logical, if TRUE then the one standard error rule is applied for selecting the optimal parameters. The one standard error rule selects the most parsimonious model having evaluation criterion no more than one standard error worse than that of the best evaluation criterion (default is FALSE).
#' @param cure_cutoff numeric value representing the cutoff time value that represents subjects not experiencing the event by this time are cured. This value is used to produce a proxy for the unobserved cure status when calculating C-statistic and AUC (default is 5 representing 5 years). Users should be careful to note the time scale of their data and adjust this according to the time scale and clinical application.
#' @param parallel logical. If TRUE, parallel processing is performed for K-fold CV using \code{foreach} and the \pkg{doMC} package is required.
#' @param seed optional integer representing the random seed. Setting the random seed fosters reproducibility of the results.
#' @param verbose logical, if TRUE running information is printed to the console (default is FALSE).
#' @param ... additional arguments.
#'
#' @return \item{b0 }{Estimated intercept for the incidence portion of the model.}
#' @return \item{b }{Estimated coefficients for the incidence portion of the model.}
#' @return \item{beta }{Estimated coefficients for the latency portion of the model.}
#' @return \item{alpha }{Estimated shape parameter if the Weibull model is fit.}
#' @return \item{rate }{Estimated rate parameter if the Weibull or exponential model is fit.}
#' @return \item{logLik }{Log-likelihood value.}
#' @return \item{selected.step.inc }{Iteration step selected for the incidence portion of the model using cross-validation. NULL when fdr.control is TRUE.}
#' @return \item{selected.step.lat }{Iteration step selected for the latency portion of the model using cross-validation. NULL when fdr.control is TRUE.}
#' @return \item{max.c }{Maximum C-statistic achieved}
#' @return \item{max.auc }{Maximum AUC for cure prediction achieved; only output when \code{measure.inc="auc"}.}
#' @return \item{selected.index.inc }{Indices of selected variables for the incidence portion of the model when \code{fdr.control=TRUE}. If none selected, \code{int(0)} will be returned.}
#' @return \item{selected.index.lat }{Indices of selected variables for the latency portion of the model when \code{fdr.control=TRUE}. If none selected, \code{int(0)} will be returned.}
#' @return \item{call}{the matched call.}
#'
#' @seealso \code{\link{curegmifs}}
#'
#' @keywords models
#' @keywords regression
#'
#' @export
#'
#' @import knockoff
#' @import doMC
#' @import stats
#' @import survival
#' @import foreach
#' @import plyr
#' @importFrom methods is
#'
#' @references Fu, H., Nicolet, D., Mrozek, K., Stone, R. M., Eisfeld, A. K., Byrd, J. C., Archer, K. J. (2022) Controlled variable selection in Weibull mixture cure models for high-dimensional data. \emph{Statistics in Medicine}, \bold{41}(22), 4340--4366.
#'
#' @seealso \code{\link{curegmifs}}
#'
#' @examples
#' library(survival)
#' set.seed(123)
#' temp <- generate_cure_data(N = 100, J = 15, nTrue = 3, A = 1.8, rho = 0.2)
#' training <- temp$Training
#'
# Fit a penalized Weibull MCM using GMIFS selecting parameters using 2-fold CV
#' fit.cv <- cv_curegmifs(Surv(Time, Censor) ~ ., data = training,
#'                       x.latency = training, fdr.control = FALSE,
#'                       maxit = 500, epsilon = 0.01,
#'                       n_folds = 2, seed = 23, verbose = TRUE)
cv_curegmifs <- function(formula, data, subset, x.latency=NULL, model="weibull",
         penalty.factor.inc=NULL, penalty.factor.lat=NULL, fdr.control=FALSE, fdr=0.2,
         epsilon=0.001, thresh=1e-5, scale=TRUE, maxit=1e4,inits=NULL,
         n_folds=5, measure.inc="c", one.se=FALSE, cure_cutoff=5,
         parallel=FALSE, seed=NULL, verbose=TRUE, ...) {
  mf <- match.call(expand.dots = FALSE)
  cl <- match.call()
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  if (m[1] == 0) stop("A \"formula\" argument is required")
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  if (missing(data))
    mf[["data"]] <- environment(formula)
  if (missing(data))
    data <- environment(formula)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  model <- c("weibull", "exponential")[pmatch(model, c("weibull", "exponential"))]
  measure.inc <- c("c","auc")[pmatch(measure.inc, c("c","auc"))]
  y <- model.response(mf)
  event <- y[, 2]
  time <- y[, 1]
  x <- model.matrix(mt, mf)
  if (!is.null(x.latency)) {
    if (missing(subset))
      r <- TRUE
    else {
      e <- substitute(subset)
      r <- eval(e, data)
      if (!is.logical(r))
        stop("'subset' must evaluate to logical")
      r <- r & !is.na(r)
    }
    if ("character" %in% is(x.latency) || "numeric" %in% is(x.latency)) {
      nl <- as.list(1:ncol(data))
      names(nl) <- names(data)
      vars <- eval(substitute(x.latency), nl, parent.frame())
      x.latency <- data[r, vars, drop = FALSE]
      x.latency <- as.matrix(x.latency)
    }
    else if ("matrix" %in% is(x.latency) || "data.frame" %in% is(x.latency)) {
      text <- parse(text=cl)[[2]]
      survnames <- strsplit(as.character(text), ",")
      time.name <- substr(survnames[[2]][1], 6, nchar(survnames[[2]][1]))
      censor.name <- trimws(strsplit(survnames[[2]][2], ")")[[1]][1])
      x.latency <- x.latency[r, !(colnames(x.latency) %in% c(time.name, censor.name)) , drop = FALSE]
      x.latency <- as.matrix(x.latency)
    }  else if ("formula" %in% is(x.latency)) {
      x.latency <- model.matrix(update.formula(x.latency, new = ~ .-1), data)
    }
  }
  x.inc <- x
  is.intercept <- grep("Intercept", dimnames(x.inc)[[2]])
  if (length(is.intercept) == 1) {
    x.inc <- x.inc[, -is.intercept, drop = FALSE]
  }
  x.lat <- x.latency
  if (nrow(x.inc) != nrow(x.lat) | nrow(x.lat) != length(time) | length(time)!= length(event))
    stop("Input dimension mismatch")
  if(class(x.inc)[1] == "data.frame" | class(x.lat)[1] == "data.frame"){
    x.inc = as.matrix(x.inc); x.lat = as.matrix(x.lat)
  }
  if (is.na(match(model, c("weibull","exponential")))) {
    stop("Error: Only 'weibull' or 'exponential' available for model parameter.")
  }
  if (is.null(penalty.factor.inc))
    penalty.factor.inc <-rep(1, ncol(x.inc))
  if (is.null(penalty.factor.lat))
    penalty.factor.lat <-rep(1, ncol(x.lat))
  if (any(!c(penalty.factor.inc, penalty.factor.inc)%in%c(0,1)))
    stop("Penalty factors specified in penalty.factor.inc and penalty.factor.inc can only include 0 or 1")
  if (fdr>1 | fdr<0)
    stop("FDR should be between 0 and 1")
  if (!is.null(inits))
    inits = inits_check(model, N=length(time), penalty.factor.inc, penalty.factor.lat, inits)
  if (is.na(match(measure.inc, c("c","auc"))))
    stop("Only 'c', 'auc' available for 'measure.inc' parameter")
  X_u = self_scale(x.inc[,penalty.factor.inc==0, drop=F], scale)
  X_p = self_scale(x.inc[,penalty.factor.inc==1, drop=F], scale)
  W_u = self_scale(x.lat[,penalty.factor.lat==0, drop=F], scale)
  W_p = self_scale(x.lat[,penalty.factor.lat==1, drop=F], scale)
  if (fdr.control) {
    res = cv.gmifs.fdr(X_u, X_p, W_u, W_p, time, event, model, fdr, thresh, nIter=maxit, epsilon, inits,
                       n_folds, measure.inc, one.se, cure_cutoff, parallel, seed, verbose)
    # KJA 05-015 added b0 = res$b0, b = res$b, beta = res$beta, rate = res$rate, alpha = res$alpha,
    b = rep(NA, ncol(x.inc))
    beta = rep(NA, ncol(x.lat))
    b[penalty.factor.inc==0] = res$b_u
    b[penalty.factor.inc==1] = res$b_p
    beta[penalty.factor.lat==0] = res$beta_u
    beta[penalty.factor.lat==1] = res$beta_p
    names(b) <- colnames(x.inc)
    names(beta) <- colnames(x.latency)
    output = list(b0 = res$b0, b = b, beta = beta, rate = res$rate, alpha = res$alpha, selected.index.inc = (1:ncol(x.inc))[penalty.factor.inc==1][res$selected_b],
                  selected.index.lat = (1:ncol(x.lat))[penalty.factor.lat==1][res$selected_beta])
    if(!is.null(colnames(x.inc)))
      names(output$selected.index.inc) <- colnames(x.inc)[output$selected.index.inc]
    if(!is.null(colnames(x.lat)))
      names(output$selected.index.lat) <- colnames(x.lat)[output$selected.index.lat]
  } else {
    res = cv.gmifs.nofdr(X_u, X_p, W_u, W_p, time, event, model, thresh, nIter=maxit, epsilon, inits,
                         n_folds, measure.inc, one.se, cure_cutoff, parallel, seed, verbose)
    b = rep(NA, ncol(x.inc))
    beta = rep(NA, ncol(x.lat))
    b[penalty.factor.inc==0] = res$b_u
    b[penalty.factor.inc==1] = res$b_p
    beta[penalty.factor.lat==0] = res$beta_u
    beta[penalty.factor.lat==1] = res$beta_p
    names(b) <- colnames(x.inc)
    names(beta) <- colnames(x.latency)
    output = list(selected.step.inc=res$model.select.inc, selected.step.lat=res$model.select.lat,
                  b0=res$b0, b=b, beta=beta, rate=res$rate,
                  logLik = res$logLik, max.c=res$max.c, x.incidence = x.inc,
                  x.latency=x.lat, y = y, model = model, scale = scale, method="GMIFS")
    if(model=="weibull") output$alpha = res$alpha
    class(output) <- "cvmixturecure"
    if(measure.inc=="auc") output$max.auc=res$max.auc
  }
  output$method = "GMIFS"
  output$model = model
  output$cv = TRUE
  output$y = y
  output$x.incidence = x.inc
  output$x.latency = x.latency
  output$scale = scale
  output$call = cl
  output$fdr.control = fdr.control
  class(output) <- "mixturecure"
  return(output)
}
