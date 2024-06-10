#' Fit penalized mixture cure model using the E-M algorithm with cross-validation for parameter tuning
#'
#' @description
#' Fits a penalized parametric and semi-parametric mixture cure model (MCM) using the E-M algorithm with with k-fold cross-validation for parameter tuning. The lasso (L1), MCP and SCAD penalty are supported for the Cox MCM while only lasso is currently supported for parametric MCMs. When FDR controlled variable selection is used, the model-X knockoffs method is applied and indices of selected variables are returned.
#'
#' @param formula an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The response must be a survival object as returned by the \code{Surv} function while the variables on the right side of the formula are the covariates that are included in the incidence portion of the model.
#' @param data a data.frame in which to interpret the variables named in the \code{formula} or in the \code{subset} argument.
#' @param subset an optional expression indicating which subset of observations to be used in the fitting process, either a numeric or factor variable should be used in subset, not a character variable. All observations are included by default.
#' @param x.latency specifies the variables to be included in the latency portion of the model and can be either a matrix of predictors, a model formula with the right hand side specifying the latency variables, or the same data.frame passed to the \code{data} parameter. Note that when using the model formula syntax for \code{x.latency} it cannot handle \code{x.latency = ~ .}.
#' @param model type of regression model to use for the latency portion of mixture cure model. Can be "cox", "weibull", or "exponential" (default is "cox").
#' @param penalty type of penalty function. Can be "lasso", "MCP", or "SCAD" (default is "lasso").
#' @param penalty.factor.inc vector of binary indicators representing the penalty to apply to each incidence coefficient: 0 implies no shrinkage and 1 implies shrinkage. If not supplied, 1 is applied to all incidence variables.
#' @param penalty.factor.lat vector of binary indicators representing the penalty to apply to each latency coefficient: 0 implies no shrinkage and 1 implies shrinkage. If not supplied, 1 is applied to all latency variables.
#' @param fdr.control logical, if TRUE, model-X knockoffs are used for FDR-controlled variable selection and indices of selected variables are returned (default is FALSE).
#' @param fdr numeric value in (0, 1) range specifying the target FDR level to use for variable selection when \code{fdr.control=TRUE} (default is 0.2).
#' @param grid.tuning logical, if TRUE a 2-D grid tuning approach is used to select the optimal pair of \eqn{\lambda_b} and \eqn{\lambda_{\beta}} penalty parameters for the incidence and latency portions of the model, respectively. Otherwise the \eqn{\lambda_b} and \eqn{\lambda_{\beta}} are selected from a 1-D sequence and are equal to one another (default is FALSE).
#' @param thresh small numeric value. The iterative process stops when the differences between successive expected penalized complete-data log-likelihoods for both incidence and latency components are less than this specified level of tolerance (default is 10^-3).
#' @param scale logical, if TRUE the predictors are centered and scaled.
#' @param maxit maximum number of passes over the data for each lambda. If not specified, 100 is applied when \code{penalty = "lasso"} and 1000 is applied when \code{penalty = "MCP"} or \code{penalty = "SCAD"}.
#' @param inits an optional list specifiying the initial value for the incidence intercept (\code{itct}), a numeric vector for the unpenalized incidence coefficients (\code{b_u}), and a numeric vector for unpenalized latency coefficients (\code{beta_u}).  For parametric models, it should also include a numeric value for the rate parameter (\code{lambda}) when \code{model = "weibull"} or \code{model = "exponential"}, and a numeric value for the shape parameter (\code{alpha}) when \code{model = "weibull"}. When \code{model = "cox"}, it should also include a numeric vector for the latency survival probabilities \eqn{S_u(t_i|w_i)} for i=1,...,N (\code{survprob}). Penalized coefficients are initialized to zero. If \code{inits} is not specified or improperly specified, initialization is automatically provided by the function.
#' @param lambda.inc.list a numeric vector used to search for the optimal \eqn{\lambda_b} tuning parameter. If not supplied, the function computes a \eqn{\lambda_b} sequence based on \code{nlambda.inc} and \code{lambda.min.ratio.inc}. If \code{grid.tuning=FALSE}, the same sequence should be used for both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}.
#' @param lambda.lat.list a numeric vector used to search for the optimal \eqn{\lambda_{\beta}} tuning parameter. If not supplied, the function computes a \eqn{\lambda_{\beta}} sequence based on \code{nlambda.lat} and \code{lambda.min.ratio.lat}. If \code{grid.tuning=FALSE}, the same sequence should be used for both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}.
#' @param nlambda.inc an integer specifying the number of values to search for the optimal \eqn{\lambda_b} tuning parameter; default is 10 if \code{grid.tuning=TRUE} and 50 otherwise.
#' @param nlambda.lat an integer specifying the number of values to search for the optimal \eqn{\lambda_{\beta}} tuning parameter; default is 10 if \code{grid.tuning=TRUE} and 50 otherwise.
#' @param gamma.inc numeric value for the penalization parameter \eqn{\gamma} for variables in the incidence portion of the model when \code{penalty = "MCP"} or \code{penalty = "SCAD"} (default is 3).
#' @param gamma.lat numeric value for the penalization parameter \eqn{\gamma} for variables in the latency portion of the model when \code{penalty = "MCP"} or \code{penalty = "SCAD"} (default is 3).
#' @param lambda.min.ratio.inc numeric value in (0,1) representing the smallest value for \eqn{\lambda_b} as a fraction of \code{lambda.max.inc}, the data-derived entry value at which essentially all penalized variables in the incidence portion of the model have a coefficient estimate of 0 (default is 0.1).
#' @param lambda.min.ratio.lat numeric value in (0.1) representing the smallest value for \eqn{\lambda_{\beta}} as a fraction of \code{lambda.max.lat}, the data-derived entry value at essentially all penalized variables in the latency portion of the model have a coefficient estimate of 0 (default is 0.1).
#' @param n_folds an integer specifying the number of folds for the k-fold cross-valiation procedure (default is 5).
#' @param measure.inc character string specifying the evaluation criterion used in selecting the optimal \eqn{\lambda_b}. Can be "c" or "auc"; default is "c". If \code{measure.inc="c"}, the C-statistic using the cure status weighting (CSW) method proposed by Asano and Hirakawa (2017) is used to select both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}. If \code{measure.inc="auc"}, the AUC for cure prediction using the mean score imputation (MSI) method proposed by Asano et al. (2014) is used to select \eqn{\lambda_b} while the C-statistic with CSW is used for \eqn{\lambda_{\beta}}.
#' @param one.se logical, if TRUE then the one standard error rule is applied for selecting the optimal parameters. The one standard error rule selects the most parsimonious model having evaluation criterion no more than one standard error worse than that of the best evaluation criterion (default is FALSE).
#' @param cure_cutoff numeric value representing the cutoff time value that represents subjects not experiencing the event by this time are cured. This value is used to produce a proxy for the unobserved cure status when calculating C-statistic and AUC (default is 5 representing 5 years). Users should be careful to note the time scale of their data and adjust this according to the time scale and clinical application.
#' @param parallel logical. If TRUE, parallel processing is performed for K-fold CV using \code{foreach} and the \pkg{doMC} package is required.
#' @param seed optional integer representing the random seed. Setting the random seed fosters reproducibility of the results.
#' @param verbose logical, if TRUE running information is printed to the console (default is FALSE).
#' @param ... additional arguments.
#'
#' @return \item{b0}{Estimated intercept for the incidence portion of the model.}
#' @return \item{b}{Estimated coefficients for the incidence portion of the model.}
#' @return \item{beta}{Estimated coefficients for the latency portion of the model.}
#' @return \item{alpha}{Estimated shape parameter if the Weibull model is fit.}
#' @return \item{rate}{Estimated rate parameter if the Weibull or exponential model is fit.}
#' @return \item{logLik.inc }{Expected penalized complete-data log-likelihood for the incidence portion of the model.}
#' @return \item{logLik.lat }{Expected penalized complete-data log-likelihood for the latency portion of the model.}
#' @return \item{selected.lambda.inc }{Value of \eqn{\lambda_b} selected using cross-validation. NULL when fdr.control is TRUE.}
#' @return \item{selected.lambda.lat }{Value of \eqn{\lambda_{\beta}} selected using cross-validation. NULL when fdr.control is TRUE.}
#' @return \item{max.c}{Maximum C-statistic achieved.}
#' @return \item{max.auc}{Maximum AUC for cure prediction achieved; only output when \code{measure.inc="auc"}.}
#' @return \item{selected.index.inc }{Indices of selected variables for the incidence portion of the model when \code{fdr.control=TRUE}. If no variables are selected, \code{int(0)} will be returned.}
#' @return \item{selected.index.lat }{Indices of selected variables for the latency portion of the model when \code{fdr.control=TRUE}. If no variables are selected, \code{int(0)} will be returned.}
#' @return \item{call}{the matched call.}

#' @export
#'
#' @import knockoff
#' @import doMC
#' @import stats
#' @import survival
#' @import glmnet
#' @importFrom methods is
#'
#' @references Archer, K. J., Fu, H., Mrozek, K., Nicolet, D., Mims, A. S., Uy, G. L., Stock, W., Byrd, J. C., Hiddemann, W., Braess, J., Spiekermann, K., Metzeler, K. H., Herold, T., Eisfeld, A.-K. (2024) Identifying long-term survivors and those at higher or lower risk of relapse among patients with cytogenetically normal acute myeloid leukemia using a high-dimensional mixture cure model. \emph{Journal of Hematology & Oncology}, \bold{17}:28.
#'
#' @seealso \code{\link{cureem}}
#'
#' @keywords models
#' @keywords regression
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(N = 200, J = 25, nTrue = 5, A = 1.8)
#' training <- temp$Training

# Fit a penalized Cox MCM selecting parameters using 2-fold CV
#' fit.cv <- cv_cureem(Surv(Time, Censor) ~ ., data = training,
#'                  x.latency = training, fdr.control = FALSE,
#'                  grid.tuning = FALSE, nlambda.inc = 10, nlambda.lat = 10,
#'                  n_folds = 2, seed = 23, verbose = TRUE)
# Select variables from a penalized Weibull MCM with FDR control and CV
#' fit.cv.fdr <- cv_cureem(Surv(Time, Censor) ~ ., data = training,
#'                  x.latency = training, model = "weibull", penalty = "lasso",
#'                  fdr.control = TRUE, grid.tuning = FALSE, nlambda.inc = 10,
#'                  nlambda.lat = 10, n_folds = 2, seed = 23, verbose = TRUE)
cv_cureem <- function(formula, data, subset, x.latency=NULL, model="cox", penalty="lasso",
         penalty.factor.inc=NULL, penalty.factor.lat=NULL,
         fdr.control=FALSE, fdr=0.2, grid.tuning = FALSE,
         thresh=1e-03, scale=TRUE, maxit=NULL, inits=NULL,
         lambda.inc.list=NULL, lambda.lat.list=NULL,
         nlambda.inc = NULL, nlambda.lat = NULL, gamma.inc=3, gamma.lat=3,
         lambda.min.ratio.inc = 0.1, lambda.min.ratio.lat = 0.1,
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
  model <- c("cox", "weibull", "exponential")[pmatch(model, c("cox", "weibull", "exponential"))]
  penalty <- c("lasso","MCP","SCAD")[pmatch(penalty, c("lasso","MCP","SCAD"))]
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
      x.latency <- model.matrix(update.formula(x.latency, new = ~ .-1, data = data), data)
    }
  }
  x.inc <- x
  is.intercept <- grep("Intercept", dimnames(x.inc)[[2]])
  if (length(is.intercept) == 1) {
    x.inc <- x.inc[, -is.intercept, drop = FALSE]
  }
  x.lat <- x.latency
  if (is.null(maxit)) maxit=ifelse(penalty=="lasso", 100, 1000)
  if (is.na(match(model, c("cox", "weibull","exponential")))) {
    stop("Error: Only 'cox', 'weibull', or 'exponential' available for model parameter.")
  }
  if (grid.tuning) {
    if (is.null(nlambda.inc)) nlambda.inc <- 10
    if (is.null(nlambda.lat)) nlambda.lat <- 10
  } else {
    if (is.null(nlambda.inc)) nlambda.inc <- 50
    if (is.null(nlambda.inc)) nlambda.lat <- 50
  }
  if (nrow(x.inc) != nrow(x.lat) | nrow(x.lat) != length(time) | length(time)!= length(event))
    stop("Input dimension mismatch")
  if (is.null(penalty.factor.inc))
    penalty.factor.inc <-rep(1, ncol(x.inc))
  if (is.null(penalty.factor.lat))
    penalty.factor.lat <-rep(1, ncol(x.lat))
  if (any(!c(penalty.factor.inc, penalty.factor.inc)%in%c(0,1)))
    stop("Penalty factors specified in penalty.factor.inc and penalty.factor.inc can only include 0 or 1")
  if (is.na(match(penalty, c("lasso","MCP","SCAD"))))
    stop("Only 'lasso', 'MCP', 'SCAD' available for 'penalty' parameter")
  if(any(c(lambda.inc.list, lambda.lat.list, gamma.inc, gamma.lat)<=0))
    stop("Penalty pamameters lambda and gamma should be positive")
  if (fdr>1 | fdr<0)
    stop("FDR should be between 0 and 1")
  if (lambda.min.ratio.inc>1 | lambda.min.ratio.inc<0 | lambda.min.ratio.lat>1 | lambda.min.ratio.lat<0)
    stop("lambda.min.ratio.inc and lambda.min.ratio.lat should be between 0 and 1")
  if (any(c(lambda.inc.list, lambda.lat.list, gamma.inc, gamma.lat)<=0))
    stop("Penalty pamameters lambda and gamma should be positive")
  if (!is.null(inits))
    inits = inits_check(model, N=length(time), penalty.factor.inc, penalty.factor.lat, inits)
  if (!measure.inc%in%c("c","auc"))
    stop("Only 'c', 'auc' available for 'measure.inc' parameter")
  if (model != "cox" & penalty != "lasso"){
    warning("MCP/SCAD penalized parametric models are not currently supported. An L1 penalized model was fitted instead.")
    penalty = "lasso"
  }
  X_u = self_scale(x.inc[,penalty.factor.inc==0, drop=F], scale)
  X_p = self_scale(x.inc[,penalty.factor.inc==1, drop=F], scale)
  W_u = self_scale(x.lat[,penalty.factor.lat==0, drop=F], scale)
  W_p = self_scale(x.lat[,penalty.factor.lat==1, drop=F], scale)
  if (fdr.control) {
    res = cv.em.fdr(X_u, X_p, W_u, W_p, time, event, model, penalty, fdr, thresh, nIter=maxit,
                    penalty.factor.inc, penalty.factor.lat,
                    grid.tuning, lambda.inc.list, lambda.lat.list, nlambda.inc, nlambda.lat,
                    lambda.min.ratio.inc, lambda.min.ratio.lat, gamma.inc, gamma.lat, inits,
                    n_folds, measure.inc, one.se, cure_cutoff, parallel, seed, verbose)
    if(!is.null(X_u)) {
      b = rep(NA, ncol(x.inc))
      b[penalty.factor.inc==0] = res$b[1:ncol(X_u)]
      b[penalty.factor.inc==1] = res$b[(ncol(X_u)+1):(length(res$b))]
    } else b = res$b
    if (!is.null(W_u)) {
      beta = rep(NA, ncol(x.lat))
      beta[penalty.factor.lat==0] = res$beta[1:ncol(W_u)]
      beta[penalty.factor.lat==1] = res$beta[(ncol(W_u)+1):(length(res$beta))]
    } else beta = res$beta
    names(b) <- colnames(x.inc)
    names(beta) <- colnames(x.latency)
    output = list(b0 = res$b0, b = b, beta = beta, rate = res$rate, alpha = res$alpha, selected.index.inc = (1:ncol(x.inc))[penalty.factor.inc==1][res$selected_b],
                  selected.index.lat = (1:ncol(x.lat))[penalty.factor.lat==1][res$selected_beta])
    if(!is.null(colnames(x.inc)))
      names(output$selected.index.inc) <- colnames(x.inc)[output$selected.index.inc]
    if(!is.null(colnames(x.lat)))
      names(output$selected.index.lat) <- colnames(x.lat)[output$selected.index.lat]
  } else {
    res = cv.em.nofdr(X_u, X_p, W_u, W_p, time, event, model, penalty, thresh, nIter=maxit,
                      grid.tuning, lambda.inc.list, lambda.lat.list, nlambda.inc, nlambda.lat,
                      lambda.min.ratio.inc, lambda.min.ratio.lat, gamma.inc, gamma.lat, inits,
                      n_folds, measure.inc, one.se, cure_cutoff, parallel, seed, verbose)
    if(!is.null(X_u)) {
      b = rep(NA, ncol(x.inc))
      b[penalty.factor.inc==0] = res$b[1:ncol(X_u)]
      b[penalty.factor.inc==1] = res$b[(ncol(X_u)+1):(length(res$b))]
    } else b = res$b
    if (!is.null(W_u)) {
      beta = rep(NA, ncol(x.lat))
      beta[penalty.factor.lat==0] = res$beta[1:ncol(W_u)]
      beta[penalty.factor.lat==1] = res$beta[(ncol(W_u)+1):(length(res$beta))]
    } else beta = res$beta
    names(b) <- colnames(x.inc)
    names(beta) <- colnames(x.latency)
    output = res
    output$b = b
    output$beta = beta
  }
  output$method = "EM"
  output$model = model
  output$penalty = penalty
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
