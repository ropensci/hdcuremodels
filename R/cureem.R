#' Fit penalized mixture cure model using the E-M algorithm
#'
#' @description
#' Fits a penalized parametric and semi-parametric mixture cure model (MCM)
#' using the E-M algorithm with user-specified penalty parameters. The lasso
#' (L1), MCP, and SCAD penalty is supported for the Cox MCM while only lasso is
#' currently supported for parametric MCMs.
#'
#' @param formula an object of class "\code{formula}" (or one that can be
#' coerced to that class): a symbolic description of the model to be fitted.
#' The response must be a survival object as returned by the \code{Surv}
#' function while the variables on the right side of the formula are the
#' covariates that are included in the incidence portion of the model.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula} or in the \code{subset} argument.
#' @param subset an optional expression indicating which subset of observations
#' to be used in the fitting process, either a numeric or factor variable
#' should be used in subset, not a character variable. All observations are
#' included by default.
#' @param x_latency specifies the variables to be included in the latency
#' portion of the model and can be either a matrix of predictors, a model
#' formula with the right hand side specifying the latency variables, or the
#' same data.frame passed to the \code{data} parameter. Note that when using
#' the model formula syntax for \code{x_latency} it cannot handle
#' \code{x_latency = ~ .}.
#' @param model type of regression model to use for the latency portion of
#' mixture cure model. Can be "cox", "weibull", or "exponential" (default is
#' "cox").
#' @param penalty type of penalty function. Can be "lasso", "MCP", or "SCAD"
#' (default is "lasso").
#' @param penalty_factor_inc vector of binary indicators representing the
#' penalty to apply to each incidence coefficient: 0 implies no shrinkage
#' and 1 implies shrinkage. If not supplied, 1 is applied to all incidence
#' variables.
#' @param penalty_factor_lat vector of binary indicators representing the
#' penalty to apply to each latency coefficient: 0 implies no shrinkage and 1
#' implies shrinkage. If not supplied, 1 is applied to all latency variables.
#' @param thresh small numeric value. The iterative process stops when the
#' differences between successive expected penalized complete-data
#' log-likelihoods for both incidence and latency components are less than this
#' specified level of tolerance (default is 10^-3).
#' @param scale logical, if TRUE the predictors are centered and scaled.
#' @param maxit integer specifying the maximum number of passes over the data
#' for each lambda. If not specified, 100 is applied when
#' \code{penalty = "lasso"} and 1000 is applied when \code{penalty = "MCP"} or
#' \code{penalty = "SCAD"}.
#' @param inits an optional list specifiying the initial value for the
#' incidence intercept (\code{itct}), a numeric vector for the unpenalized
#' incidence coefficients (\code{b_u}), and a numeric vector for unpenalized
#' latency coefficients (\code{beta_u}).  For parametric models, it should also
#' include a numeric value for the rate parameter (\code{lambda}) when
#' \code{model = "weibull"} or \code{model = "exponential"}, and a numeric
#' value for the shape parameter (\code{alpha}) when \code{model = "weibull"}.
#' When \code{model = "cox"}, it should also include a numeric vector for the
#' latency survival probabilities \eqn{S_u(t_i|w_i)} for i=1,...,N
#' (\code{survprob}). Penalized coefficients are initialized to zero.
#' If \code{inits} is not specified or improperly specified, initialization is
#' automatically provided by the function.
#' @param lambda_inc numeric value for the penalization parameter \eqn{\lambda}
#' for variables in the incidence portion of the model.
#' @param lambda_lat numeric value for the penalization parameter \eqn{\lambda}
#' for variables in the latency portion of the model.
#' @param gamma_inc numeric value for the penalization parameter \eqn{\gamma}
#' for variables in the incidence portion of the model when
#' \code{penalty = "MCP"} or \code{penalty = "SCAD"} (default is 3).
#' @param gamma_lat numeric value for the penalization parameter \eqn{\gamma}
#' for variables in the latency portion of the model when \code{penalty = "MCP"}
#' or \code{penalty = "SCAD"} (default is 3).
#' @param ... additional arguments.
#'
#' @return \item{b_path}{Matrix representing the solution path of the
#' coefficients in the incidence portion of the model. Row is step and column
#' is variable.}
#' @return \item{beta_path}{Matrix representing the solution path of lthe
#' coefficients in the latency portion of the model. Row is step and column
#' is variable.}
#' @return \item{b0_path }{Vector representing the solution path of the
#' intercept in the incidence portion of the model.}
#' @return \item{logLik_inc }{Vector representing the expected penalized
#' complete-data log-likelihood for the incidence portion of the model for
#' each step in the solution path.}
#' @return \item{logLik_lat }{Vector representing the expected penalized
#' complete-data log-likelihood for the latency portion of the model for each
#' step in the solution path.}
#' @return \item{x_incidence}{Matrix representing the design matrix of the
#' incidence predictors.}
#' @return \item{x_latency}{Matrix representing the design matrix of the
#' latency predictors.}
#' @return \item{y}{Vector representing the survival object response as
#' returned by the \code{Surv} function }
#' @return \item{model}{Character string indicating the type of regression
#' model used for the latency portion of mixture cure model ("weibull" or
#' "exponential").}
#' @return \item{scale}{Logical value indicating whether the predictors were
#' centered and scaled.}
#' @return \item{method}{Character string indicating the EM alogoritm was used
#' in fitting the mixture cure model.}
#' @return \item{rate_path}{Vector representing the solution path of the rate
#' parameter for the Weibull or exponential density in the latency portion of
#' the model.}
#' @return \item{alpha_path}{Vector representing the solution path of the shape
#' parameter for the Weibull density in the latency portion of the model.}
#' @return \item{call}{the matched call.}
#' @export
#'
#' @import stats
#' @import survival
#' @import glmnet
#' @import foreach
#' @import plyr
#' @importFrom methods is
#'
#' @references Archer, K. J., Fu, H., Mrozek, K., Nicolet, D., Mims, A. S.,
#' Uy, G. L., Stock, W., Byrd, J. C., Hiddemann, W., Braess, J.,
#' Spiekermann, K., Metzeler, K. H., Herold, T., Eisfeld, A.-K. (2024)
#' Identifying long-term survivors and those at higher or lower risk of relapse
#' among patients with cytogenetically normal acute myeloid leukemia using a
#' high-dimensional mixture cure model. \emph{Journal of Hematology & Oncology},
#' \bold{17}:28.
#'
#' @srrstats {G1.0} *Statistical Software should list at least one primary reference from published academic literature.*
#' @srrstats {G1.1} *Statistical Software should document whether the algorithm(s) it implements are: An improvement on other implementations of similar algorithms in **R**.
#'
#' @seealso \code{\link{cv_cureem}}
#'
#' @keywords models
#' @keywords regression
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 80, j = 100, n_true = 10, a = 1.8)
#' training <- temp$training
#' fit <- cureem(Surv(Time, Censor) ~ .,
#'   data = training, x_latency = training,
#'   model = "cox", penalty = "lasso", lambda_inc = 0.1,
#'   lambda_lat = 0.1, gamma_inc = 6, gamma_lat = 10
#' )
cureem <- function(formula, data, subset, x_latency = NULL,
                   model = c("cox", "weibull", "exponential"),
                   penalty = c("lasso", "MCP", "SCAD"),
                   penalty_factor_inc = NULL, penalty_factor_lat = NULL,
                   thresh = 1e-03, scale = TRUE, maxit = NULL, inits = NULL,
                   lambda_inc = 0.1, lambda_lat = 0.1, gamma_inc = 3,
                   gamma_lat = 3, ...) {
  mf <- match.call(expand.dots = FALSE)
  cl <- match.call()
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  if (m[1] == 0) stop("A \"formula\" argument is required")
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  if (missing(data)) {
    mf[["data"]] <- environment(formula)
  }
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  model <- match.arg(model)
  penalty <- match.arg(penalty)
  y <- model.response(mf)
  event <- y[, 2]
  time <- y[, 1]
  x <- model.matrix(mt, mf)
  if (!is.null(x_latency)) {
    if (missing(subset)) {
      r <- TRUE
    } else {
      e <- substitute(subset)
      r <- eval(e, data)
      if (!is.logical(r)) {
        stop("'subset' must evaluate to logical")
      }
      r <- r & !is.na(r)
    }
    if ("character" %in% is(x_latency) || "numeric" %in% is(x_latency)) {
      nl <- as.list(seq_len(ncol(data)))
      names(nl) <- names(data)
      vars <- eval(substitute(x_latency), nl, parent.frame())
      x_latency <- data[r, vars, drop = FALSE]
      x_latency <- as.matrix(x_latency)
    } else if ("matrix" %in% is(x_latency) || "data.frame" %in% is(x_latency)) {
      text <- parse(text = cl)[[2]]
      survnames <- strsplit(as.character(text), ",")
      time_name <- substr(survnames[[2]][1], 6, nchar(survnames[[2]][1]))
      censor_name <- trimws(strsplit(survnames[[2]][2], ")")[[1]][1])
      x_latency <- x_latency[r, !(colnames(x_latency) %in% c(
        time_name,
        censor_name
      )),
      drop = FALSE
      ]
      x_latency <- as.matrix(x_latency)
    } else if ("formula" %in% is(x_latency)) {
      x_latency <- model.matrix(update.formula(x_latency, new = ~ . - 1), data)
    }
  }
  x_inc <- x
  is_intercept <- grep("Intercept", dimnames(x_inc)[[2]])
  if (length(is_intercept) == 1) {
    x_inc <- x_inc[, -is_intercept, drop = FALSE]
  }
  x_lat <- x_latency
  if (is.null(maxit)) maxit <- ifelse(penalty == "lasso", 100, 1000)
  if (nrow(x_inc) != nrow(x_lat) || nrow(x_lat) != length(time) || length(time) != length(event)) {
    stop("Input dimension mismatch")
  }
  if (class(x_inc)[1] == "data.frame" || class(x_lat)[1] == "data.frame") {
    x_inc <- as.matrix(x_inc)
    x_lat <- as.matrix(x_lat)
  }
  if (is.null(penalty_factor_inc)) {
    penalty_factor_inc <- rep(1, ncol(x_inc))
  }
  if (is.null(penalty_factor_lat)) {
    penalty_factor_lat <- rep(1, ncol(x_lat))
  }
  if (any(!c(penalty_factor_inc, penalty_factor_inc) %in% c(0, 1))) {
    stop("Penalty factors specified in penalty_factor_inc and
         penalty_factor_inc can only include 0 or 1")
  }
  if (any(c(lambda_inc, lambda_lat, gamma_inc, gamma_lat) <= 0)) {
    stop("Penalty pamameters lambda and gamma should be positive")
  }
  if (!is.null(inits)) {
    inits <- inits_check(model,
      N = length(time), penalty_factor_inc,
      penalty_factor_lat, inits
    )
  }
  if (model != "cox" && penalty != "lasso") {
    warning("MCP/SCAD penalized parametric models are not currently supported.
            An L1 penalized model was fitted instead.")
    penalty <- "lasso"
  }
  x_u <- self_scale(x_inc[, penalty_factor_inc == 0, drop = FALSE], scale)
  x_p <- self_scale(x_inc[, penalty_factor_inc == 1, drop = FALSE], scale)
  w_u <- self_scale(x_lat[, penalty_factor_lat == 0, drop = FALSE], scale)
  w_p <- self_scale(x_lat[, penalty_factor_lat == 1, drop = FALSE], scale)
  if (model == "cox" && penalty == "lasso") {
    fit <- cox_l1(
      x_u, x_p, w_u, w_p, time, event, lambda_inc, lambda_lat,
      inits, maxit, thresh
    )
  } else if (model == "cox" && penalty %in% c("MCP", "SCAD")) {
    fit <- cox_mcp_scad(
      x_u, x_p, w_u, w_p, time, event, penalty, lambda_inc,
      lambda_lat, gamma_inc, gamma_lat, inits, maxit, thresh
    )
  } else if (model == "weibull") {
    fit <- weib_EM(
      x_u, x_p, w_u, w_p, time, event, lambda_inc, lambda_lat,
      inits, maxit, thresh
    )
  } else if (model == "exponential") {
    fit <- exp_EM(
      x_u, x_p, w_u, w_p, time, event, lambda_inc, lambda_lat,
      inits, maxit, thresh
    )
  }
  b_path <- matrix(NA, nrow = dim(fit$b_p_path)[1], ncol = dim(x_inc)[2])
  beta_path <- matrix(NA, nrow = dim(fit$beta_p_path)[1], ncol = dim(x_lat)[2])
  b_path[, penalty_factor_inc == 0] <- fit$b_u_path
  b_path[, penalty_factor_inc == 1] <- fit$b_p_path
  beta_path[, penalty_factor_lat == 0] <- fit$beta_u_path
  beta_path[, penalty_factor_lat == 1] <- fit$beta_p_path
  colnames(b_path) <- colnames(x_inc)
  colnames(beta_path) <- colnames(x_lat)
  b0_path <- fit$itct_path
  logLik_inc <- fit$lik_inc
  logLik_lat <- fit$lik_lat
  output <- list(
    b_path = b_path, beta_path = beta_path,
    b0_path = b0_path, logLik_inc = logLik_inc,
    logLik_lat = logLik_lat, x_incidence = x_inc,
    x_latency = x_lat, y = y, model = model, scale = scale,
    method = "EM", call = cl
  )
  if (model %in% c("exponential", "weibull")) {
    output$rate <- fit$lambda_path
  }
  if (model == "weibull") {
    output$alpha <- fit$alpha_path
  }
  output$cv <- FALSE
  class(output) <- "mixturecure"
  return(output)
}
