#' Fit penalized parametric mixture cure model using the GMIFS algorithm
#'
#' @description
#' Fits a penalized Weibull or exponential mixture cure model using the generalized monotone incremental forward stagewise (GMIFS) algorithm and yields solution paths for parameters in the incidence and latency portions of the model.
#'
#' @param formula an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The response must be a survival object as returned by the \code{Surv} function while the variables on the right side of the formula are the covariates that are included in the incidence portion of the model.
#' @param data a data.frame in which to interpret the variables named in the \code{formula} or in the \code{subset} argument.
#' @param subset an optional expression indicating which subset of observations to be used in the fitting process, either a numeric or factor variable should be used in subset, not a character variable. All observations are included by default.
#' @param x.latency specifies the variables to be included in the latency portion of the model and can be either a matrix of predictors, a model formula with the right hand side specifying the latency variables, or the same data.frame passed to the \code{data} parameter. Note that when using the model formula syntax for \code{x.latency} it cannot handle \code{x.latency = ~ .}.
#' @param model type of regression model to use for the latency portion of mixture cure model. Can be "weibull" or "exponential"; default is "weibull".
#' @param penalty.factor.inc vector of binary indicators representing the penalty to apply to each incidence coefficient: 0 implies no shrinkage and 1 implies shrinkage. If not supplied, 1 is applied to all incidence variables.
#' @param penalty.factor.lat vector of binary indicators representing the penalty to apply to each latency coefficient: 0 implies no shrinkage and 1 implies shrinkage. If not supplied, 1 is applied to all latency variables.
#' @param epsilon small numeric value reflecting the incremental value used to update a coefficient at a given step (default is 0.001).
#' @param thresh small numeric value. The iterative process stops when the differences between successive expected penalized complete-data log-likelihoods for both incidence and latency components are less than this specified level of tolerance (default is 10^-5).
#' @param scale logical, if TRUE the predictors are centered and scaled.
#' @param maxit integer specifying the maximum number of steps to run in the iterative algorithm (default is 10^4).
#' @param inits an optional list specifiying the initial value for the incidence intercept (\code{itct}), a numeric vector for the unpenalized incidence coefficients (\code{b_u}), and a numeric vector for unpenalized latency coefficients (\code{beta_u}), a numeric value for the rate parameter (\code{lambda}), and a numeric value for the shape parameter (\code{alpha}) when \code{model = "weibull"}. If not supplied or improperly supplied, initialization is automatically provided by the function.
#' @param verbose logical, if TRUE running information is printed to the console (default is FALSE).
#' @param ... additional arguments.
#'
#' @return \item{b_path}{Matrix representing the solution path of the coefficients in the incidence portion of the model. Row is step and column is variable.}
#' @return \item{beta_path}{Matrix representing the solution path of lthe coefficients in the latency portion of the model. Row is step and column is variable.}
#' @return \item{b0_path}{Vector representing the solution path of the intercept in the incidence portion of the model.}
#' @return \item{rate_path}{Vector representing the solution path of the rate parameter for the Weibull or exponential density in the latency portion of the model.}
#' @return \item{logLik}{Vector representing the log-likelihood for each step in the solution path.}
#' @return \item{x.incidence}{Matrix representing the design matrix of the incidence predictors.}
#' @return \item{x.latency}{Matrix representing the design matrix of the latency predictors.}
#' @return \item{y}{Vector representing the survival object response as returned by the \code{Surv} function }
#' @return \item{model}{Character string indicating the type of regression model used for the latency portion of mixture cure model ("weibull" or "exponential").}
#' @return \item{scale}{Logical value indicating whether the predictors were centered and scaled.}
#' @return \item{alpha_path}{Vector representing the solution path of the shape parameter for the Weibull density in the latency portion of the model.}
#' @return \item{call}{the matched call.}
#'
#' @references Fu, H., Nicolet, D., Mrozek, K., Stone, R. M., Eisfeld, A. K., Byrd, J. C., Archer, K. J. (2022) Controlled variable selection in Weibull mixture cure models for high-dimensional data. \emph{Statistics in Medicine}, \bold{41}(22), 4340--4366.
#'
#' @seealso \code{\link{cv_curegmifs}}
#' @export
#'
#' @import stats
#' @import survival
#' @importFrom methods is
#'
#' @keywords models
#' @keywords regression
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(N = 100, J = 10, nTrue = 10, A = 1.8)
#' training <- temp$Training
#'
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'          data = training, x.latency = training,
#'          model = "weibull", thresh = 1e-4, maxit = 2000, epsilon = 0.01,
#'          verbose = FALSE)
curegmifs <- function(formula, data, subset, x.latency=NULL, model="weibull", penalty.factor.inc=NULL, penalty.factor.lat=NULL,
         epsilon=0.001, thresh=1e-5, scale=TRUE, maxit=1e4,inits=NULL, verbose=TRUE, ...) {
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
  if(any(!c(penalty.factor.inc, penalty.factor.inc)%in%c(0,1)))
    stop("Penalty factors specified in penalty.factor.inc and penalty.factor.inc can only include 0 or 1")
  if(!is.null(inits))
    inits = inits_check(model, N=length(time), penalty.factor.inc, penalty.factor.lat, inits)
  X_u = self_scale(x.inc[,penalty.factor.inc == 0, drop = F], scale)
  X_p = self_scale(x.inc[,penalty.factor.inc == 1, drop = F], scale)
  W_u = self_scale(x.lat[,penalty.factor.lat == 0, drop = F], scale)
  W_p = self_scale(x.lat[,penalty.factor.lat == 1, drop = F], scale)
  if(model=="exponential")
    res = exp_cure(X_u, X_p, W_u, W_p, time, event, epsilon, thresh, maxit, inits, verbose)
  else if(model=="weibull")
    res = weibull.cure(X_u, X_p, W_u, W_p, time, event, epsilon, thresh, maxit, inits, verbose)
  nstep = length(res$logLikelihood)
  if(nstep==maxit) warning("Maximum step of iterations achieved. Algorithm may not converge.")
  b_path = matrix(NA, nstep, ncol(x.inc))
  beta_path = matrix(NA, nstep, ncol(x.lat))
  b_path[,penalty.factor.inc==0] = res$b_u_path
  b_path[,penalty.factor.inc==1] = res$b_p_path
  beta_path[,penalty.factor.lat==0] = res$beta_u_path
  beta_path[,penalty.factor.lat==1] = res$beta_p_path
  colnames(b_path) <- colnames(x.inc)
  colnames(beta_path) <- colnames(x.lat)
  output = list(b_path = b_path, beta_path = beta_path, b0_path = res$itct_path,
                rate_path = res$lambda_path, logLik = res$logLikelihood, x.incidence = x.inc,
                x.latency=x.lat, y = y, model = model, scale = scale, method="GMIFS", call = cl)
  if(model=="weibull") output$alpha_path = res$alpha_path
  output$cv = FALSE
  class(output) <- "mixturecure"
  return(output)
}
