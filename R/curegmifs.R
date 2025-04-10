#' Fit penalized parametric mixture cure model using the GMIFS algorithm
#'
#' @description
#' Fits a penalized Weibull or exponential mixture cure model using the
#' generalized monotone incremental forward stagewise (GMIFS) algorithm and
#' yields solution paths for parameters in the incidence and latency portions of
#' the model.
#'
#' @param formula an object of class "\code{formula}" (or one that can be
#' coerced to that class): a symbolic description of the model to be fitted. The
#' response must be a survival object as returned by the \code{Surv} function
#' while the variables on the right side of the formula are the covariates that
#' are included in the incidence portion of the model.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula} or in the \code{subset} argument.
#' @param subset an optional expression indicating which subset of observations
#' to be used in the fitting process, either a numeric or factor variable should
#' be used in subset, not a character variable. All observations are included by
#' default.
#' @param x_latency specifies the variables to be included in the latency
#' portion of the model and can be either a matrix of predictors, a model
#' formula with the right hand side specifying the latency variables, or the
#' same data.frame passed to the \code{data} parameter. Note that when using the
#' model formula syntax for \code{x_latency} it cannot handle
#' \code{x_latency = ~ .}.
#' @param model type of regression model to use for the latency portion of
#' mixture cure model. Can be "weibull" or "exponential"; default is "weibull".
#' @param penalty_factor_inc vector of binary indicators representing the
#' penalty to apply to each incidence coefficient: 0 implies no shrinkage and 1
#' implies shrinkage. If not supplied, 1 is applied to all incidence variables.
#' @param penalty_factor_lat vector of binary indicators representing the
#' penalty to apply to each latency coefficient: 0 implies no shrinkage and 1
#' implies shrinkage. If not supplied, 1 is applied to all latency variables.
#' @param epsilon small numeric value reflecting the incremental value used to
#' update a coefficient at a given step (default is 0.001).
#' @param thresh small numeric value. The iterative process stops when the
#' differences between successive expected penalized
#' log-likelihoods for both incidence and latency components are less than this
#' specified level of tolerance (default is 10^-5).
#' @param scale logical, if TRUE the predictors are centered and scaled.
#' @param maxit integer specifying the maximum number of steps to run in the
#' iterative algorithm (default is 10^4).
#' @param inits an optional list specifying the initial value for the incidence
#' intercept (\code{itct}), a numeric vector for the unpenalized incidence
#' coefficients (\code{b_u}), and a numeric vector for unpenalized latency
#' coefficients (\code{beta_u}), a numeric value for the rate parameter
#' (\code{lambda}), and a numeric value for the shape parameter (\code{alpha})
#' when \code{model = "weibull"}. If not supplied or improperly supplied,
#' initialization is automatically provided by the function.
#' @param verbose logical, if TRUE running information is printed to the
#' console (default is FALSE).
#' @param suppress_warning logical, if TRUE, suppresses echoing the warning that
#' the maximum number of iterations was reached so that the algorithm may not
#' have converged. Instead, warning is returned as part of the output with this
#' message.
#' @param ... additional arguments.
#'
#' @return \item{b_path}{Matrix representing the solution path of the
#' coefficients in the incidence portion of the model. Row is step and
#' column is variable.}
#' @return \item{beta_path}{Matrix representing the solution path of the
#' coefficients in the latency portion of the model. Row is step and column is
#' variable.}
#' @return \item{b0_path}{Vector representing the solution path of the intercept
#' in the incidence portion of the model.}
#' @return \item{rate_path}{Vector representing the solution path of the rate
#' parameter for the Weibull or exponential density in the latency portion of
#' the model.}
#' @return \item{logLik}{Vector representing the log-likelihood for each step
#' in the solution path.}
#' @return \item{x_incidence}{Matrix representing the design matrix of the
#' incidence predictors.}
#' @return \item{x_latency}{Matrix representing the design matrix of the latency
#' predictors.}
#' @return \item{y}{Vector representing the survival object response as returned
#' by the \code{Surv} function }
#' @return \item{model}{Character string indicating the type of regression model
#' used for the latency portion of mixture cure model ("weibull" or
#' "exponential").}
#' @return \item{scale}{Logical value indicating whether the predictors were
#' centered and scaled.}
#' @return \item{alpha_path}{Vector representing the solution path of the shape
#' parameter for the Weibull density in the latency portion of the model.}
#' @return \item{call}{the matched call.}
#' @return \item{warning}{message indicating whether the maximum number of
#' iterations was achieved which may indicate the model did not converge.}
#'
#' @references Fu, H., Nicolet, D., Mrozek, K., Stone, R. M., Eisfeld, A. K.,
#' Byrd, J. C., Archer, K. J. (2022) Controlled variable selection in Weibull
#' mixture cure models for high-dimensional data. \emph{Statistics in Medicine},
#' \bold{41}(22), 4340--4366.
#'
#' @srrstats {G1.0} *Statistical Software should list at least one primary reference from published academic literature.*
#' @srrstats {G1.1} *Statistical Software should document whether the algorithm(s) it implements are: The first implementation of a novel algorithm *
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {G1.3} *All statistical terminology should be clarified and unambiguously defined.*
#' @srrstats {G1.5} *Software should include all code necessary to reproduce results which form the basis of performance claims made in associated publications.*
#' @srrstats {G2.0a} *Provide explicit secondary documentation of any expectations on lengths of inputs*
#' @srrstats {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @srrstats {G2.1a} *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
#' @srrstats {G2.2} *Appropriately prohibit or restrict submission of multivariate input to parameters expected to be univariate.*
#' @srrstats {G2.3} *For univariate character input:*
#' @srrstats {G2.3b} *Either: use `tolower()` or equivalent to ensure input of character parameters is not case dependent; or explicitly document that parameters are strictly case-sensitive.*
#' @srrstats {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstats {G2.4e} *explicit conversion from factor via `as...()` functions*
#' @srrstats {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.*
#' @srrstats {G5.2} *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
#' @srrstats {G5.2a} *Every message produced within R code by `stop()`, `warning()`, `message()`, or equivalent should be unique*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstats {RE1.0} *Regression Software should enable models to be specified via a formula interface, unless reasons for not doing so are explicitly documented.*
#' @srrstats {RE1.1} *Regression Software should document how formula interfaces are converted to matrix representations of input data.*
#' @srrstats {RE1.2} *Regression Software should document expected format (types or classes) for inputting predictor variables, including descriptions of types or classes which are not accepted.*
#' @srrstats {RE1.3} *Regression Software which passes or otherwise transforms aspects of input data onto output structures should ensure that those output structures retain all relevant aspects of input data, notably including row and column names, and potentially information from other `attributes()`.*
#' @srrstats {RE1.3a} *Where otherwise relevant information is not transferred, this should be explicitly documented.*
#' @srrstats {RE1.4} *Regression Software should document any assumptions made with regard to input data; for example distributional assumptions, or assumptions that predictor data have mean values of zero. Implications of violations of these assumptions should be both documented and tested.*
#' @srrstats {RE2.3} *Where applicable, Regression Software should enable data to be centred (for example, through converting to zero-mean equivalent values; or to z-scores) or offset (for example, to zero-intercept equivalent values) via additional parameters, with the effects of any such parameters clearly documented and tested.*
#' @srrstats {RE3.0} *Issue appropriate warnings or other diagnostic messages for models which fail to converge.*
#' @srrstats {RE3.2} *Ensure that convergence thresholds have sensible default values, demonstrated through explicit documentation.*
#' @srrstats {RE3.3} *Allow explicit setting of convergence thresholds, unless reasons against doing so are explicitly documented.*
#' @srrstats {RE4.0} *Regression Software should return some form of "model" object, generally through using or modifying existing class structures for model objects (such as `lm`, `glm`, or model objects from other packages), or creating a new class of model objects.*
#' @srrstats {RE4.4} *The specification of the model, generally as a formula (via `formula()`)*
#' @srrstats {RE4.7} *Where appropriate, convergence statistics*
#' @srrstats {RE4.8} *Response variables, and associated "metadata" where applicable.*
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
#' withr::local_seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#'
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'   data = training, x_latency = training,
#'   model = "weibull", thresh = 1e-4, maxit = 2000, epsilon = 0.01,
#'   verbose = FALSE
#' )
curegmifs <- function(formula, data, subset, x_latency = NULL,
                      model = c("weibull", "exponential"),
                      penalty_factor_inc = NULL, penalty_factor_lat = NULL,
                      epsilon = 0.001, thresh = 1e-5, scale = TRUE, maxit = 1e4,
                      inits = NULL, verbose = TRUE, suppress_warning = FALSE,
                      ...) {
  mf <- match.call(expand.dots = FALSE)
  cl <- match.call()
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  if (m[1] == 0) stop("Error: A \"formula\" argument is required")
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  if (missing(data)) {
    mf[["data"]] <- environment(formula)
  }
  if (missing(data)) {
    data <- environment(formula)
  }
  if (thresh < 0) {
    stop("Error: \"thresh\" must be non-negative")
  }
  if (epsilon < 0) {
    stop("Error: \"epsilon\" must be non-negative")
  }
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  model <- tolower(model)
  model <- match.arg(model)
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
        stop("Error: 'subset' must evaluate to logical")
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
      x_latency <- x_latency[r, !(colnames(x_latency) %in%
        c(time_name, censor_name)), drop = FALSE]
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
  if (nrow(x_inc) != nrow(x_lat) || nrow(x_lat) != length(time) ||
    length(time) != length(event)) {
    stop("Error: Input dimension mismatch")
  }
  if (class(x_inc)[1] == "data.frame" || class(x_lat)[1] == "data.frame") {
    x_inc <- as.matrix(x_inc)
    x_lat <- as.matrix(x_lat)
  }
  if (is.na(match(model, c("weibull", "exponential")))) {
    stop("Error: Only 'weibull' or 'exponential' available for model
         parameter.")
  }
  if (is.null(penalty_factor_inc)) {
    penalty_factor_inc <- rep(1, ncol(x_inc))
  }
  if (is.null(penalty_factor_lat)) {
    penalty_factor_lat <- rep(1, ncol(x_lat))
  }
  if (any(!c(class(penalty_factor_inc), class(penalty_factor_lat)) %in% "numeric")) {
    stop("Error: Penalty factors specified in penalty_factor_inc and penalty_factor_lat
         must be numeric vectors comprised of 0 or 1")
  }
  if (any(!c(penalty_factor_inc, penalty_factor_lat) %in% c(0, 1))) {
    stop("Error: Penalty factors specified in penalty_factor_inc and penalty_factor_lat
         can only include 0 or 1")
  }
  if (!is.null(inits)) {
    inits <- inits_check(model,
      N = length(time), penalty_factor_inc,
      penalty_factor_lat, inits
    )
  }
  x_u <- self_scale(x_inc[, penalty_factor_inc == 0, drop = FALSE], scale)
  x_p <- self_scale(x_inc[, penalty_factor_inc == 1, drop = FALSE], scale)
  w_u <- self_scale(x_lat[, penalty_factor_lat == 0, drop = FALSE], scale)
  w_p <- self_scale(x_lat[, penalty_factor_lat == 1, drop = FALSE], scale)
  if (model == "exponential") {
    res <- exp_cure(
      x_u, x_p, w_u, w_p, time, event, epsilon, thresh, maxit,
      inits, verbose
    )
  } else if (model == "weibull") {
    res <- weibull.cure(
      x_u, x_p, w_u, w_p, time, event, epsilon, thresh, maxit,
      inits, verbose
    )
  }
  nstep <- length(res$logLikelihood)
  warning <- "No warnings"
  if (nstep == maxit) {
    warning <-
      "Warning: Maximum step of iterations achieved.
                           Algorithm may not converge."
    if (!suppress_warning)
      warning("Warning: Maximum step of iterations achieved.
                           Algorithm may not converge.")
  }
  b_path <- matrix(NA, nrow = nstep, ncol = ncol(x_inc))
  beta_path <- matrix(NA, nrow = nstep, ncol = ncol(x_lat))
  b_path[, penalty_factor_inc == 0] <- res$b_u_path
  b_path[, penalty_factor_inc == 1] <- res$b_p_path
  beta_path[, penalty_factor_lat == 0] <- res$beta_u_path
  beta_path[, penalty_factor_lat == 1] <- res$beta_p_path
  colnames(b_path) <- colnames(x_inc)
  colnames(beta_path) <- colnames(x_lat)
  output <- list(
    b_path = b_path, beta_path = beta_path,
    b0_path = res$itct_path, rate_path = res$lambda_path,
    logLik = res$logLikelihood, x_incidence = x_inc,
    x_latency = x_lat, y = y, model = model, scale = scale,
    method = "GMIFS", call = cl
  )
  if (model == "weibull") output$alpha_path <- res$alpha_path
  output$cv <- FALSE
  output$warning <- warning
  class(output) <- "mixturecure"
  return(output)
}
