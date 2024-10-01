#' Fit penalized mixture cure model using the E-M algorithm with
#' cross-validation for parameter tuning
#'
#' @description
#' Fits a penalized parametric and semi-parametric mixture cure model (MCM)
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
#' the \code{formula} or in the \code{subset} argument.
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
#' use for variable selection when \code{fdr_control=TRUE} (default is 0.2).
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
#' @param inits an optional list specifiying the initial value for the incidence
#' intercept (\code{itct}), a numeric vector for the unpenalized incidence
#' coefficients (\code{b_u}), and a numeric vector for unpenalized latency
#' coefficients (\code{beta_u}).  For parametric models, it should also include
#' a numeric value for the rate parameter (\code{lambda}) when
#' \code{model = "weibull"} or \code{model = "exponential"}, and a numeric
#' value for the shape parameter (\code{alpha}) when \code{model = "weibull"}.
#' When \code{model = "cox"}, it should also include a numeric vector for the
#' latency survival probabilities \eqn{S_u(t_i|w_i)} for i=1,...,N
#' (\code{survprob}). Penalized coefficients are initialized to zero. If
#' \code{inits} is not specified or improperly specified, initialization is
#' automatically provided by the function.
#' @param lambda_inc_list a numeric vector used to search for the optimal
#' \eqn{\lambda_b} tuning parameter. If not supplied, the function computes a
#' \eqn{\lambda_b} sequence based on \code{nlambda_inc} and
#' \code{lambda_min_ratio_inc}. If \code{grid_tuning=FALSE}, the same sequence
#' should be used for both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}.
#' @param lambda_lat_list a numeric vector used to search for the optimal
#' \eqn{\lambda_{\beta}} tuning parameter. If not supplied, the function
#' computes a \eqn{\lambda_{\beta}} sequence based on \code{nlambda_lat} and
#' \code{lambda_min_ratio_lat}. If \code{grid_tuning=FALSE}, the same sequence
#' should be used for both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}.
#' @param nlambda_inc an integer specifying the number of values to search for
#' the optimal \eqn{\lambda_b} tuning parameter; default is 10 if
#' \code{grid_tuning=TRUE} and 50 otherwise.
#' @param nlambda_lat an integer specifying the number of values to search
#' for the optimal \eqn{\lambda_{\beta}} tuning parameter; default is 10 if
#' \code{grid_tuning=TRUE} and 50 otherwise.
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
#' in selecting the optimal \eqn{\lambda_b}. Can be "c" or "auc"; default is
#' "c". If \code{measure_inc="c"}, the C-statistic using the cure status
#' weighting (CSW) method proposed by Asano and Hirakawa (2017) is used to
#' select both \eqn{\lambda_b} and \eqn{\lambda_{\beta}}. If
#' \code{measure_inc="auc"}, the AUC for cure prediction using the mean score
#' imputation (MSI) method proposed by Asano et al. (2014) is used to select
#' \eqn{\lambda_b} while the C-statistic with CSW is used for
#' \eqn{\lambda_{\beta}}.
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

#' @seealso \code{\link{cureem}}
#'
#' @keywords models
#' @keywords regression
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 200, j = 25, n_true = 5, a = 1.8)
#' training <- temp$training

# Fit a penalized Cox MCM selecting parameters using 2-fold CV
#' fit.cv <- cv_cureem(Surv(Time, Censor) ~ ., data = training,
#'                  x_latency = training, fdr_control = FALSE,
#'                  grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
#'                  n_folds = 2, seed = 23, verbose = TRUE)
# Select variables from a penalized Weibull MCM with FDR control and CV
#' fit.cv.fdr <- cv_cureem(Surv(Time, Censor) ~ ., data = training,
#'                  x_latency = training, model = "weibull", penalty = "lasso",
#'                  fdr_control = TRUE, grid_tuning = FALSE, nlambda_inc = 10,
#'                  nlambda_lat = 10, n_folds = 2, seed = 23, verbose = TRUE)
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
                      verbose = TRUE, ...) {
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
  model <- match.arg(tolower(model))
  penalty <- match.arg(tolower(penalty))
  measure_inc <- match.arg(tolower(measure_inc))
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
      if (!is.logical(r))
        stop("'subset' must evaluate to logical")
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
    }  else if ("formula" %in% is(x_latency)) {
      x_latency <- model.matrix(update.formula(x_latency, new = ~ . - 1,
                                               data = data), data)
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
      length(time) != length(event))
    stop("Input dimension mismatch")
  if (is.null(penalty_factor_inc))
    penalty_factor_inc <- rep(1, ncol(x_inc))
  if (is.null(penalty_factor_lat))
    penalty_factor_lat <- rep(1, ncol(x_lat))
  if (any(!c(penalty_factor_inc, penalty_factor_inc) %in% c(0, 1)))
    stop("Penalty factors specified in penalty_factor_inc and penalty_factor_inc
         can only include 0 or 1")
  if (any(c(lambda_inc_list, lambda_lat_list, gamma_inc, gamma_lat) <= 0))
    stop("Penalty pamameters lambda and gamma should be positive")
  if (fdr > 1 || fdr < 0)
    stop("FDR should be between 0 and 1")
  if (lambda_min_ratio_inc > 1 || lambda_min_ratio_inc < 0 ||
        lambda_min_ratio_lat > 1 || lambda_min_ratio_lat < 0)
    stop("lambda_min_ratio_inc and lambda_min_ratio_lat should be between
         0 and 1")
  if (any(c(lambda_inc_list, lambda_lat_list, gamma_inc, gamma_lat) <= 0))
    stop("Penalty pamameters lambda and gamma should be positive")
  if (!is.null(inits))
    inits <- inits_check(model, N = length(time), penalty_factor_inc,
                         penalty_factor_lat, inits)
  if (model != "cox" && penalty != "lasso") {
    warning("MCP/SCAD penalized parametric models are not currently supported.
            An L1 penalized model was fitted instead.")
    penalty <- "lasso"
  }
  x_u <- self_scale(x_inc[, penalty_factor_inc == 0, drop = FALSE], scale)
  x_p <- self_scale(x_inc[, penalty_factor_inc == 1, drop = FALSE], scale)
  w_u <- self_scale(x_lat[, penalty_factor_lat == 0, drop = FALSE], scale)
  w_p <- self_scale(x_lat[, penalty_factor_lat == 1, drop = FALSE], scale)
  if (fdr_control) {
    res <- cv.em.fdr(x_u, x_p, w_u, w_p, time, event, model, penalty, fdr,
                     thresh, nIter = maxit, penalty_factor_inc,
                     penalty_factor_lat, grid_tuning, lambda_inc_list,
                     lambda_lat_list, nlambda_inc, nlambda_lat,
                     lambda_min_ratio_inc, lambda_min_ratio_lat,
                     gamma_inc, gamma_lat, inits, n_folds, measure_inc, one_se,
                     cure_cutoff, parallel, seed, verbose)
    if (!is.null(x_u)) {
      b <- rep(NA, ncol(x_inc))
      b[penalty_factor_inc == 0] <- res$b[seq_len(ncol(x_u))]
      b[penalty_factor_inc == 1] <- res$b[(ncol(x_u) + 1):(length(res$b))]
    } else {
      b <- res$b
    }
    if (!is.null(w_u)) {
      beta <- rep(NA, ncol(x_lat))
      beta[penalty_factor_lat == 0] <- res$beta[seq_len(ncol(w_u))]
      beta[penalty_factor_lat == 1] <- res$beta[(ncol(w_u) + 1):
                                                  (length(res$beta))]
    } else {
      beta <- res$beta
    }
    names(b) <- colnames(x_inc)
    names(beta) <- colnames(x_latency)
    output <- list(b0 = res$b0, b = b, beta = beta, rate = res$rate,
                   alpha = res$alpha,
                   selected_index_inc =
                     (seq_len(ncol(x_inc)))[penalty_factor_inc == 1][
                                                                res$selected_b],
                   selected_index_lat =
                     (seq_len(ncol(x_lat)))[penalty_factor_lat == 1][
                                                             res$selected_beta])
    if (!is.null(colnames(x_inc)))
      names(output$selected_index_inc) <-
        colnames(x_inc)[output$selected_index_inc]
    if (!is.null(colnames(x_lat)))
      names(output$selected_index_lat) <-
        colnames(x_lat)[output$selected_index_lat]
  } else {
    res <- cv.em.nofdr(x_u, x_p, w_u, w_p, time, event, model, penalty,
                       thresh, nIter = maxit, grid_tuning, lambda_inc_list,
                       lambda_lat_list, nlambda_inc, nlambda_lat,
                       lambda_min_ratio_inc, lambda_min_ratio_lat, gamma_inc,
                       gamma_lat, inits, n_folds, measure_inc, one_se,
                       cure_cutoff, parallel, seed, verbose)
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
