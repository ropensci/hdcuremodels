#' Fit a penalized parametric mixture cure model using the GMIFS algorithm with
#' cross-validation for model selection
#'
#' @description
#' Fits a penalized Weibull or exponential mixture cure model using the
#' generalized monotone incremental forward stagewise (GMIFS) algorithm with
#' k-fold cross-validation to select the optimal iteration step along the
#' solution path. When FDR controlled variable selection is used, the model-X
#' knockoffs method is applied and indices of selected variables are returned.
#'
#' @param formula an object of class "\code{formula}" (or one that can be
#' coerced to that class): a symbolic description of the model to be fitted.
#' The response must be a survival object as returned by the \code{Surv}
#' function while the variables on the right side of the formula are the
#' covariates that are included in the incidence portion of the model.
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
#' @param fdr_control logical, if TRUE, model-X knockoffs are used for
#' FDR-controlled variable selection and indices of selected variables are
#' returned (default is FALSE).
#' @param fdr numeric value in (0, 1) range specifying the target FDR level to
#' use for variable selection when \code{fdr_control=TRUE} (default is 0.2).
#' @param epsilon small numeric value reflecting incremental value used to
#' update a coefficient at a given step (default is 0.001).
#' @param thresh small numeric value. The iterative process stops when the
#' differences between successive expected penalized complete-data
#' log-likelihoods for both incidence and latency components are less than this
#' specified level of tolerance (default is 10^-5).
#' @param scale logical, if TRUE the predictors are centered and scaled.
#' @param maxit integer specifying the maximum number of steps to run in the
#' iterative algorithm (default is 10^4).
#' @param inits an optional list specifiying the initial value for the incidence
#' intercept (\code{itct}), a numeric vector for the unpenalized incidence
#' coefficients (\code{b_u}), and a numeric vector for unpenalized latency
#' coefficients (\code{beta_u}), a numeric value for the rate parameter
#' (\code{lambda}), and a numeric value for the shape parameter (\code{alpha})
#' when \code{model = "weibull"}. If not supplied or improperly supplied,
#' initialization is automatically provided by the function.
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
#' @return \item{b0 }{Estimated intercept for the incidence portion of the
#' model.}
#' @return \item{b }{Estimated coefficients for the incidence portion of the
#' model.}
#' @return \item{beta }{Estimated coefficients for the latency portion of the
#' model.}
#' @return \item{alpha }{Estimated shape parameter if the Weibull model is fit.}
#' @return \item{rate }{Estimated rate parameter if the Weibull or exponential
#' model is fit.}
#' @return \item{logLik }{Log-likelihood value.}
#' @return \item{selected.step.inc }{Iteration step selected for the incidence
#' portion of the model using cross-validation. NULL when fdr_control is TRUE.}
#' @return \item{selected.step.lat }{Iteration step selected for the latency
#' portion of the model using cross-validation. NULL when fdr_control is TRUE.}
#' @return \item{max.c }{Maximum C-statistic achieved}
#' @return \item{max.auc }{Maximum AUC for cure prediction achieved; only output
#' when \code{measure_inc="auc"}.}
#' @return \item{selected_index_inc }{Indices of selected variables for the
#' incidence portion of the model when \code{fdr_control=TRUE}. If none
#' selected, \code{int(0)} will be returned.}
#' @return \item{selected_index_lat }{Indices of selected variables for the
#' latency portion of the model when \code{fdr_control=TRUE}. If none selected,
#' \code{int(0)} will be returned.}
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
#' @import doParallel
#' @import parallel
#' @import stats
#' @import survival
#' @import foreach
#' @import plyr
#' @importFrom methods is
#'
#' @references Fu, H., Nicolet, D., Mrozek, K., Stone, R. M., Eisfeld, A. K.,
#' Byrd, J. C., Archer, K. J. (2022) Controlled variable selection in Weibull
#' mixture cure models for high-dimensional data. \emph{Statistics in Medicine},
#' \bold{41}(22), 4340--4366.
#'
#' @seealso \code{\link{curegmifs}}
#'
#' @examples
#' library(survival)
#' set.seed(123)
#' temp <- generate_cure_data(n = 100, j = 15, n_true = 3, a = 1.8, rho = 0.2)
#' training <- temp$training
#'
# Fit a penalized Weibull MCM using GMIFS selecting parameters using 2-fold CV
#' fit.cv <- cv_curegmifs(Surv(Time, Censor) ~ ., data = training,
#'                       x_latency = training, fdr_control = FALSE,
#'                       maxit = 450, epsilon = 0.01 ,n_folds = 2,
#'                      seed = 23, verbose = TRUE)
cv_curegmifs <- function(formula, data, subset, x_latency = NULL,
                         model = c("weibull", "exponential"),
                         penalty_factor_inc = NULL, penalty_factor_lat = NULL,
                         fdr_control = FALSE, fdr = 0.2, epsilon = 0.001,
                         thresh = 1e-5, scale = TRUE, maxit = 1e4, inits = NULL,
                         n_folds = 5, measure_inc = c("c", "auc"),
                         one_se = FALSE, cure_cutoff = 5, parallel = FALSE,
                         seed = NULL, verbose = TRUE, ...) {
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
  model <- match.arg(model)
  measure_inc <- match.arg(measure_inc)
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
                                    c(time_name, censor_name)) , drop = FALSE]
      x_latency <- as.matrix(x_latency)
    }  else if ("formula" %in% is(x_latency)) {
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
        length(time) != length(event))
    stop("Input dimension mismatch")
  if (class(x_inc)[1] == "data.frame" || class(x_lat)[1] == "data.frame") {
    x_inc <- as.matrix(x_inc)
    x_lat <- as.matrix(x_lat)
  }
  if (is.null(penalty_factor_inc))
    penalty_factor_inc <- rep(1, ncol(x_inc))
  if (is.null(penalty_factor_lat))
    penalty_factor_lat <- rep(1, ncol(x_lat))
  if (any(!c(penalty_factor_inc, penalty_factor_inc) %in% c(0, 1)))
    stop("Penalty factors specified in penalty_factor_inc and penalty_factor_inc
         can only include 0 or 1")
  if (fdr > 1 || fdr < 0)
    stop("FDR should be between 0 and 1")
  if (!is.null(inits))
    inits <- inits_check(model, N = length(time), penalty_factor_inc,
                         penalty_factor_lat, inits)
  x_u <- self_scale(x_inc[, penalty_factor_inc == 0, drop = FALSE], scale)
  x_p <- self_scale(x_inc[, penalty_factor_inc == 1, drop = FALSE], scale)
  w_u <- self_scale(x_lat[, penalty_factor_lat == 0, drop = FALSE], scale)
  w_p <- self_scale(x_lat[, penalty_factor_lat == 1, drop = FALSE], scale)
  if (fdr_control) {
    res <- cv.gmifs.fdr(x_u, x_p, w_u, w_p, time, event, model, fdr, thresh,
                        nIter = maxit, epsilon, inits, n_folds, measure_inc,
                        one_se, cure_cutoff, parallel, seed, verbose)
    b <- rep(NA, ncol(x_inc))
    beta <- rep(NA, ncol(x_lat))
    b[penalty_factor_inc == 0] <- res$b_u
    b[penalty_factor_inc == 1] <- res$b_p
    beta[penalty_factor_lat == 0] <- res$beta_u
    beta[penalty_factor_lat == 1] <- res$beta_p
    names(b) <- colnames(x_inc)
    names(beta) <- colnames(x_latency)
    output <- list(b0 = res$b0, b = b, beta = beta, rate = res$rate,
                   alpha = res$alpha, selected_index_inc =
                     (seq_len(ncol(x_inc)))[penalty_factor_inc
                                            == 1][res$selected_b],
                   selected_index_lat =
                     (seq_len(ncol(x_lat)))[penalty_factor_lat
                                            == 1][res$selected_beta])
    if (!is.null(colnames(x_inc)))
      names(output$selected_index_inc) <-
        colnames(x_inc)[output$selected_index_inc]
    if (!is.null(colnames(x_lat)))
      names(output$selected_index_lat) <-
        colnames(x_lat)[output$selected_index_lat]
  } else {
    res <- cv.gmifs.nofdr(x_u, x_p, w_u, w_p, time, event, model, thresh,
                          nIter = maxit, epsilon, inits, n_folds, measure_inc,
                          one_se, cure_cutoff, parallel, seed, verbose)
    b <- rep(NA, ncol(x_inc))
    beta <- rep(NA, ncol(x_lat))
    b[penalty_factor_inc == 0] <- res$b_u
    b[penalty_factor_inc == 1] <- res$b_p
    beta[penalty_factor_lat == 0] <- res$beta_u
    beta[penalty_factor_lat == 1] <- res$beta_p
    names(b) <- colnames(x_inc)
    names(beta) <- colnames(x_latency)
    output <- list(selected_step_inc = res$model.select.inc,
                   selected_step_lat = res$model.select.lat,
                   b0 = res$b0, b = b, beta  = beta, rate = res$rate,
                   logLik = res$logLik, max_c = res$max.c, x_incidence = x_inc,
                   x_latency = x_lat, y = y, model = model, scale = scale,
                   method = "GMIFS")
    if (model == "weibull") output$alpha <- res$alpha
    class(output) <- "cvmixturecure"
    if (measure_inc == "auc") output$max.auc <- res$max.auc
  }
  output$method <- "GMIFS"
  output$model <- model
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
