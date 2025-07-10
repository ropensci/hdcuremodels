#' Fit penalized mixture cure model using the E-M algorithm with
#' cross-validation for parameter tuning
#'
#' @description
#' Fits penalized parametric and semi-parametric mixture cure models (MCM)
#' using the E-M algorithm with with k-fold cross-validation for parameter
#' tuning. The lasso (L1), MCP and SCAD penalty are supported for the Cox MCM
#' while only lasso is currently supported for parametric MCMs. When FDR
#' controlled variable selection is used, the model-X knockoffs method is
#' applied and indices of selected variables are returned.
#'
#' @param formula an object of class "\code{formula}" (or one that can be
#' coerced to that class): a symbolic description of the model to be fitted.
#' The response must be a survival object as returned by the \code{Surv}
#' function while the variables on the right side of the formula are the
#' covariates that are included in the incidence portion of the model.
#' @param data a data.frame in which to interpret the variables named in
#' the \code{formula} or in the \code{subset} argument. Rows with missing data are
#' omitted (only \code{na.action = na.omit} is operational) therefore users may
#' want to impute missing data prior to calling this function.
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
#' penalty to apply to each incidence coefficient: 0 implies no shrinkage and
#' 1 implies shrinkage. If not supplied, 1 is applied to all incidence
#' variables.
#' @param penalty_factor_lat vector of binary indicators representing the
#' penalty to apply to each latency coefficient: 0 implies no shrinkage and 1
#' implies shrinkage. If not supplied, 1 is applied to all latency variables.
#' @param fdr_control logical, if TRUE, model-X knockoffs are used for
#' FDR-controlled variable selection and indices of selected variables are
#' returned (default is FALSE).
#' @param fdr numeric value in (0, 1) range specifying the target FDR level to
#' use for variable selection when \code{fdr_control = TRUE} (default is 0.2).
#' @param grid_tuning logical, if TRUE a 2-D grid tuning approach is used to
#' select the optimal pair of \eqn{\lambda_b} and \eqn{\lambda_{\beta}} penalty
#' parameters for the incidence and latency portions of the model, respectively.
#' Otherwise the \eqn{\lambda_b} and \eqn{\lambda_{\beta}} are selected from a
#' 1-D sequence and are equal to one another (default is FALSE).
#' @param thresh small numeric value. The iterative process stops when the
#' differences between successive expected penalized complete-data
#' log-likelihoods for both incidence and latency components are less than this
#' specified level of tolerance (default is 10^-3).
#' @param scale logical, if TRUE the predictors are centered and scaled.
#' @param maxit maximum number of passes over the data for each lambda. If not
#' specified, 100 is applied when \code{penalty = "lasso"} and 1000 is applied
#' when \code{penalty = "MCP"} or \code{penalty = "SCAD"}.
#' @param inits an optional list specifiying the initial values to be used for
#' model fitting as follows:
#' \itemize{
#' \item \code{itct} the incidence intercept.
#' \item \code{b_u} a numeric vector for the unpenalized.
#' incidence coefficients for the incidence portion of the model.
#' \item \code{beta_u} a numeric vector for unpenalized
#' latency coefficients in the incidence portion of the model.
#' \item \code{lambda} a numeric value for the rate parameter when fitting
#' either a Weibull or exponential MCM using \code{model = "weibull"} or
#' \code{model = "exponential"}.
#' \item \code{alpha} a numeric value for the shape parameter when fitting a
#' Weibull MCM using \code{model = "weibull"}.
#' \item \code{survprob} a numeric vector for the
#' latency survival probabilities \eqn{S_u(t_i|w_i)} for i=1,...,N when fitting
#' a Cox MCM \code{model = "cox"}.
#' }
#' Penalized coefficients are initialized to zero. If \code{inits} is not specified or improperly specified, initialization is
#' automatically provided by the function.
#' @param lambda_inc_list a numeric vector used to search for the optimal
#' \eqn{\lambda_b} tuning parameter. If not supplied, the function computes a
#' \eqn{\lambda_b} sequence based on \code{nlambda_inc} and
#' \code{lambda_min_ratio_inc}. If \code{grid_tuning = FALSE}, the same sequence
#' should be used for both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}.
#' @param lambda_lat_list a numeric vector used to search for the optimal
#' \eqn{\lambda_{\beta}} tuning parameter. If not supplied, the function
#' computes a \eqn{\lambda_{\beta}} sequence based on \code{nlambda_lat} and
#' \code{lambda_min_ratio_lat}. If \code{grid_tuning = FALSE}, the same sequence
#' should be used for both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}.
#' @param nlambda_inc an integer specifying the number of values to search for
#' the optimal \eqn{\lambda_b} tuning parameter; default is 10 if
#' \code{grid_tuning = TRUE} and 50 otherwise.
#' @param nlambda_lat an integer specifying the number of values to search
#' for the optimal \eqn{\lambda_{\beta}} tuning parameter; default is 10 if
#' \code{grid_tuning = TRUE} and 50 otherwise.
#' @param gamma_inc numeric value for the penalization parameter \eqn{\gamma}
#' for variables in the incidence portion of the model when
#' \code{penalty = "MCP"} or \code{penalty = "SCAD"} (default is 3).
#' @param gamma_lat numeric value for the penalization parameter \eqn{\gamma}
#' for variables in the latency portion of the model when \code{penalty = "MCP"}
#' or \code{penalty = "SCAD"} (default is 3).
#' @param lambda_min_ratio_inc numeric value in (0,1) representing the smallest
#' value for \eqn{\lambda_b} as a fraction of \code{lambda.max_inc}, the
#' data-derived entry value at which essentially all penalized variables in the
#' incidence portion of the model have a coefficient estimate of 0 (default is
#' 0.1).
#' @param lambda_min_ratio_lat numeric value in (0.1) representing the smallest
#' value for \eqn{\lambda_{\beta}} as a fraction of \code{lambda.max_lat}, the
#' data-derived entry value at essentially all penalized variables in the
#' latency portion of the model have a coefficient estimate of 0 (default is
#' 0.1).
#' @param n_folds an integer specifying the number of folds for the k-fold
#' cross-valiation procedure (default is 5).
#' @param measure_inc character string specifying the evaluation criterion used
#' in selecting the optimal \eqn{\lambda_b} which can be either
#' \itemize{
#' \item \code{"c"} specifying to use the C-statistic for cure status
#' weighting (CSW) method proposed by Asano and Hirakawa (2017) to
#' select both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}
#' \item \code{"auc"} specifying to use the AUC for cure prediction using the
#' mean score imputation (MSI) method proposed by Asano et al. (2014) to select
#' \eqn{\lambda_b} while the C-statistic with CSW is used for
#' \eqn{\lambda_{\beta}}.
#' }
#' @param one_se logical, if TRUE then the one standard error rule is applied
#' for selecting the optimal parameters. The one standard error rule selects the
#' most parsimonious model having evaluation criterion no more than one standard
#' error worse than that of the best evaluation criterion (default is FALSE).
#' @param cure_cutoff numeric value representing the cutoff time value that
#' represents subjects not experiencing the event by this time are cured. This
#' value is used to produce a proxy for the unobserved cure status when
#' calculating C-statistic and AUC (default is 5 representing 5 years). Users
#' should be careful to note the time scale of their data and adjust this
#' according to the time scale and clinical application.
#' @param parallel logical. If TRUE, parallel processing is performed for K-fold
#' CV using \code{foreach} and the \pkg{doParallel} package is required.
#' @param seed optional integer representing the random seed. Setting the random
#' seed fosters reproducibility of the results.
#' @param verbose logical, if TRUE running information is printed to the console
#' (default is FALSE).
#' @param na.action this function requires complete data so \code{"na.omit"} is
#' invoked. Users can impute missing data as an alternative prior to model fitting.
#' @param ... additional arguments.
#'
#' @return \item{b0}{Estimated intercept for the incidence portion of the
#' model.}
#' @return \item{b}{Estimated coefficients for the incidence portion of the
#' model.}
#' @return \item{beta}{Estimated coefficients for the latency portion of the
#' model.}
#' @return \item{alpha}{Estimated shape parameter if the Weibull model is fit.}
#' @return \item{rate}{Estimated rate parameter if the Weibull or exponential
#' model is fit.}
#' @return \item{logLik_inc}{Expected penalized complete-data log-likelihood for
#' the incidence portion of the model.}
#' @return \item{logLik_lat}{Expected penalized complete-data log-likelihood for
#' the latency portion of the model.}
#' @return \item{selected_lambda_inc}{Value of \eqn{\lambda_b} selected using
#' cross-validation. NULL when fdr_control is TRUE.}
#' @return \item{selected_lambda_lat}{Value of \eqn{\lambda_{\beta}} selected
#' using cross-validation. NULL when fdr_control is TRUE.}
#' @return \item{max_c}{Maximum C-statistic achieved.}
#' @return \item{max_auc}{Maximum AUC for cure prediction achieved; only output
#' when \code{measure_inc="auc"}.}
#' @return \item{selected_index_inc }{Indices of selected variables for the
#' incidence portion of the model when \code{fdr_control=TRUE}. If no variables
#' are selected, \code{int(0)} will be returned.}
#' @return \item{selected_index_lat }{Indices of selected variables for the
#' latency portion of the model when \code{fdr_control=TRUE}. If no variables
#' are selected, \code{int(0)} will be returned.}
#' @return \item{call}{the matched call.}

#' @export
#'
#' @import knockoff
#' @import doParallel
#' @import parallel
#' @import stats
#' @import survival
#' @import glmnet
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
#' @srrstats {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @srrstats {G2.1a} *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
#' @srrstats {G2.2} *Appropriately prohibit or restrict submission of multivariate input to parameters expected to be univariate.*
#' @srrstats {G2.3} *For univariate character input:*
#' @srrstats {G2.3b} *Either: use `tolower()` or equivalent to ensure input of character parameters is not case dependent; or explicitly document that parameters are strictly case-sensitive.*
#' @srrstats {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstats {G2.4e} *explicit conversion from factor via `as...()` functions*
#' @srrstats {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.*
#' @srrstats {G2.13} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#' @srrstats {G5.2} *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
#' @srrstats {G5.2a} *Every message produced within R code by `stop()`, `warning()`, `message()`, or equivalent should be unique*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstats {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
#' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstats {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*
#' @srrstats {RE1.0} *Regression Software should enable models to be specified via a formula interface, unless reasons for not doing so are explicitly documented.*
#' @srrstats {RE1.1} *Regression Software should document how formula interfaces are converted to matrix representations of input data.*
#' @srrstats {RE1.2} *Regression Software should document expected format (types or classes) for inputting predictor variables, including descriptions of types or classes which are not accepted.*
#' @srrstats {RE1.3} *Regression Software which passes or otherwise transforms aspects of input data onto output structures should ensure that those output structures retain all relevant aspects of input data, notably including row and column names, and potentially information from other `attributes()`.*
#' @srrstats {RE1.3a} *Where otherwise relevant information is not transferred, this should be explicitly documented.*
#' @srrstats {RE1.4} *Regression Software should document any assumptions made with regard to input data; for example distributional assumptions, or assumptions that predictor data have mean values of zero. Implications of violations of these assumptions should be both documented and tested.*
#' @srrstats {RE2.3} *Where applicable, Regression Software should enable data to be centred (for example, through converting to zero-mean equivalent values; or to z-scores) or offset (for example, to zero-intercept equivalent values) via additional parameters, with the effects of any such parameters clearly documented and tested.*
#' @srrstats {RE4.0} *Regression Software should return some form of "model" object, generally through using or modifying existing class structures for model objects (such as `lm`, `glm`, or model objects from other packages), or creating a new class of model objects.*
#' @srrstats {RE4.4} *The specification of the model, generally as a formula (via `formula()`)*
#' @srrstats {RE4.8} *Response variables, and associated "metadata" where applicable.*
#' @seealso \code{\link{cureem}}
#'
#' @keywords models
#' @keywords regression
#'
#' @examples
#' library(survival)
#' withr::local_seed(1234)
#' temp <- generate_cure_data(n = 200, j = 25, n_true = 5, a = 1.8)
#' training <- temp$training
# Fit a penalized Cox MCM selecting parameters using 2-fold CV
#' fit.cv <- cv_cureem(Surv(Time, Censor) ~ .,
#'   data = training,
#'   x_latency = training, fdr_control = FALSE,
#'   grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
#'   n_folds = 2, seed = 23, verbose = TRUE
#' )
# Select variables from a penalized Weibull MCM with FDR control and CV
#' fit.cv.fdr <- cv_cureem(Surv(Time, Censor) ~ .,
#'   data = training,
#'   x_latency = training, model = "weibull", penalty = "lasso",
#'   fdr_control = TRUE, grid_tuning = FALSE, nlambda_inc = 10,
#'   nlambda_lat = 10, n_folds = 2, seed = 23, verbose = TRUE)
cv_cureem <- function(formula, data, subset, x_latency = NULL,
                      model = c("cox", "weibull", "exponential"),
                      penalty = c("lasso", "MCP", "SCAD"),
                      penalty_factor_inc = NULL, penalty_factor_lat = NULL,
                      fdr_control = FALSE, fdr = 0.2, grid_tuning = FALSE,
                      thresh = 1e-03, scale = TRUE, maxit = NULL, inits = NULL,
                      lambda_inc_list = NULL, lambda_lat_list = NULL,
                      nlambda_inc = NULL, nlambda_lat = NULL, gamma_inc = 3,
                      gamma_lat = 3,
                      lambda_min_ratio_inc = 0.1, lambda_min_ratio_lat = 0.1,
                      n_folds = 5, measure_inc = c("c", "auc"), one_se = FALSE,
                      cure_cutoff = 5, parallel = FALSE, seed = NULL,
                      verbose = TRUE, na.action = na.omit, ...) {
  mf <- match.call(expand.dots = FALSE)
  cl <- match.call()
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  if (m[1] == 0) stop("Error: A \"formula\" argument is required")
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
  omitted <- attr(mf, "na.action")
  model <- match.arg(model)
  penalty <- match.arg(penalty)
  measure_inc <- match.arg(measure_inc)
  y <- model.response(mf)
  event <- y[, 2]
  time <- y[, 1]
  x <- model.matrix(mt, mf)
  if (!is.null(x_latency)) {
    if (missing(subset)) {
      r <- rep(TRUE, dim(data)[1])
      r[omitted] <- FALSE
    } else {
      e <- substitute(subset)
      r <- eval(e, data)
      if (!is.logical(r)) {
        stop("Error: 'subset' must evaluate to logical")
      }
      r <- r & !is.na(r)
      r[omitted] <- FALSE
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
      x_latency <- model.matrix(update.formula(x_latency,
        new = ~ . - 1,
        data = data
      ), data)
    }
  }
  x_inc <- x
  is_intercept <- grep("Intercept", dimnames(x_inc)[[2]])
  if (length(is_intercept) == 1) {
    x_inc <- x_inc[, -is_intercept, drop = FALSE]
  }
  x_lat <- x_latency
  if (is.null(maxit)) maxit <- ifelse(penalty == "lasso", 100, 1000)
  if (grid_tuning) {
    if (is.null(nlambda_inc)) nlambda_inc <- 10
    if (is.null(nlambda_lat)) nlambda_lat <- 10
  } else {
    if (is.null(nlambda_inc)) nlambda_inc <- 50
    if (is.null(nlambda_inc)) nlambda_lat <- 50
  }
  if (nrow(x_inc) != nrow(x_lat) || nrow(x_lat) != length(time) ||
    length(time) != length(event)) {
    stop("Error: Input dimension mismatch")
  }
  if (is.null(penalty_factor_inc)) {
    penalty_factor_inc <- rep(1, ncol(x_inc))
  }
  if (is.null(penalty_factor_lat)) {
    penalty_factor_lat <- rep(1, ncol(x_lat))
  }
  if (any(!c(class(penalty_factor_inc), class(penalty_factor_lat)) %in% "numeric")) {
    stop("Error: Penalty factors specified in penalty_factor_inc and penalty_factor_inc
         must be numeric vectors comprised of 0 or 1")
  }
  if (any(!c(penalty_factor_inc, penalty_factor_inc) %in% c(0, 1))) {
    stop("Error: Penalty factors specified in penalty_factor_inc and penalty_factor_inc
         can only include 0 or 1")
  }
  if (any(c(lambda_inc_list, lambda_lat_list, gamma_inc, gamma_lat) <= 0)) {
    stop("Error: Penalty pamameters lambda and gamma should be positive")
  }
  if (fdr > 1 || fdr < 0) {
    stop("Error: FDR should be between 0 and 1")
  }
  if (lambda_min_ratio_inc > 1 || lambda_min_ratio_inc < 0 ||
    lambda_min_ratio_lat > 1 || lambda_min_ratio_lat < 0) {
    stop("Error: lambda_min_ratio_inc and lambda_min_ratio_lat should be between
         0 and 1")
  }
  if (any(c(lambda_inc_list, lambda_lat_list, gamma_inc, gamma_lat) <= 0)) {
    stop("Error: Penalty pamameters lambda and gamma should be positive")
  }
  if (!is.null(inits)) {
    inits <- inits_check(model,
      N = length(time), penalty_factor_inc,
      penalty_factor_lat, inits
    )
  }
  if (model != "cox" && penalty != "lasso") {
    warning("Warning: MCP/SCAD penalized parametric models are not currently supported.
            An L1 penalized model was fitted instead.")
    penalty <- "lasso"
  }
  x_u <- self_scale(x_inc[, penalty_factor_inc == 0, drop = FALSE], scale)
  x_p <- self_scale(x_inc[, penalty_factor_inc == 1, drop = FALSE], scale)
  w_u <- self_scale(x_lat[, penalty_factor_lat == 0, drop = FALSE], scale)
  w_p <- self_scale(x_lat[, penalty_factor_lat == 1, drop = FALSE], scale)
  if (fdr_control) {
    res <- cv.em.fdr(x_u, x_p, w_u, w_p, time, event, model, penalty, fdr,
      thresh,
      nIter = maxit, penalty_factor_inc,
      penalty_factor_lat, grid_tuning, lambda_inc_list,
      lambda_lat_list, nlambda_inc, nlambda_lat,
      lambda_min_ratio_inc, lambda_min_ratio_lat,
      gamma_inc, gamma_lat, inits, n_folds, measure_inc, one_se,
      cure_cutoff, parallel, seed, verbose
    )
    if (!is.null(x_u)) {
      b <- rep(NA, ncol(x_inc))
      b[penalty_factor_inc == 0] <- res$b[1:ncol(x_u)]
      b[penalty_factor_inc == 1] <- res$b[(ncol(x_u) + 1):(length(res$b))]
    } else {
      b <- res$b
    }
    if (!is.null(w_u)) {
      beta <- rep(NA, ncol(x_lat))
      beta[penalty_factor_lat == 0] <- res$beta[1:ncol(w_u)]
      beta[penalty_factor_lat == 1] <- res$beta[(ncol(w_u) + 1):
      (length(res$beta))]
    } else {
      beta <- res$beta
    }
    names(b) <- colnames(x_inc)
    names(beta) <- colnames(x_latency)
    output <- list(
      b0 = res$b0, b = b, beta = beta, rate = res$rate,
      alpha = res$alpha,
      selected_index_inc =
        (seq_len(ncol(x_inc)))[penalty_factor_inc == 1][
          res$selected_b
        ],
      selected_index_lat =
        (seq_len(ncol(x_lat)))[penalty_factor_lat == 1][
          res$selected_beta
        ]
    )
    if (!is.null(colnames(x_inc))) {
      names(output$selected_index_inc) <-
        colnames(x_inc)[output$selected_index_inc]
    }
    if (!is.null(colnames(x_lat))) {
      names(output$selected_index_lat) <-
        colnames(x_lat)[output$selected_index_lat]
    }
  } else {
    res <- cv.em.nofdr(x_u, x_p, w_u, w_p, time, event, model, penalty,
      thresh,
      nIter = maxit, grid_tuning, lambda_inc_list,
      lambda_lat_list, nlambda_inc, nlambda_lat,
      lambda_min_ratio_inc, lambda_min_ratio_lat, gamma_inc,
      gamma_lat, inits, n_folds, measure_inc, one_se,
      cure_cutoff, parallel, seed, verbose
    )
    if (!is.null(x_u)) {
      b <- rep(NA, ncol(x_inc))
      b[penalty_factor_inc == 0] <- res$b[1:ncol(x_u)]
      b[penalty_factor_inc == 1] <- res$b[(ncol(x_u) + 1):(length(res$b))]
    } else {
      b <- res$b
    }
    if (!is.null(w_u)) {
      beta <- rep(NA, ncol(x_lat))
      beta[penalty_factor_lat == 0] <- res$beta[1:ncol(w_u)]
      beta[penalty_factor_lat == 1] <-
        res$beta[(ncol(w_u) + 1):(length(res$beta))]
    } else {
      beta <- res$beta
    }
    names(b) <- colnames(x_inc)
    names(beta) <- colnames(x_latency)
    output <- res
    output$b <- b
    output$beta <- beta
  }
  output$method <- "EM"
  output$model <- model
  output$penalty <- penalty
  output$cv <- TRUE
  output$y <- y
  output$x_incidence <- x_inc
  output$x_latency <- x_latency
  output$scale <- scale
  output$call <- cl
  output$fdr_control <- fdr_control
  class(output) <- "mixturecure"
  return(output)
}
