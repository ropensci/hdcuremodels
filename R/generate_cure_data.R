#' Simulate data under a mixture cure model
#'
#' @description
#' Simulate data under a mixture cure model.
#'
#' @param n an integer denoting the total sample size.
#' @param j an integer denoting the number of penalized predictors which is the
#' same for both the incidence and latency portions of the model.
#' @param nonp an integer denoting the number of unpenalized
#' predictors (which is the same for both the incidence and latency portions of
#' the model).
#' @param train_prop a numeric value in [0, 1) representing the fraction of `n`
#' to be used in forming the training dataset.
#' @param n_true an integer less than `j` denoting the number of variables truly
#' associated with the outcome (i.e., the number of covariates with nonzero
#' parameter values) among the penalized predictors.
#' @param a a numeric value denoting the effect size (signal amplitude) which is
#' the same for both the incidence and latency portions of the model.
#' @param rho a numeric value in [0, 1) representing the correlation between
#' adjacent covariates in the same block.
#' @param itct_mean a numeric value representing the expectation of the
#' incidence intercept which controls the cure rate.
#' @param cens_ub a numeric value representing the upper bound on the censoring
#' time distribution which follows a uniform distribution on (0, \code{cens_ub}].
#' @param alpha a numeric value representing the shape parameter in the Weibull
#' density.
#' @param lambda a numeric value representing the rate parameter in the Weibull
#' density.
#' @param same_signs logical, if TRUE the incidence and latency coefficients
#' have the same signs.
#' @param model type of regression model to use for the latency portion of
#' mixture cure model. Can be one of the following:
#' \itemize{
#' \item \code{"weibull"} to generate times from a Weibull distribution.
#' \item \code{"GG"} to generate times from a generalized gamma distribution.
#' \item \code{"Gompertz"} to generate times from a Gomertz distribution.
#' \item \code{"nonparametric"} to generate times non-parametrically.
#' \item  \code{"GG_baseline"} to generate times from a generalized gamma
#' baseline distribution.
#' }
#' @return \item{training}{training data.frame which includes Time, Censor, and
#' covariates. Variables prefixed with \code{"U"} indicates unpenalized
#' covariates and is equal to the value
#' passed to \code{nonp} (default is 2). Variables prefixed with \code{"X"}
#' indicates penalized covariates and is equal to the value passed to
#' \code{j}.
#' }
#' @return \item{training_y}{the true status for the training set: uncured = 1;
#' cured = 0}
#' @return \item{testing}{testing data.frame which includes Time, Censor,  Y
#' (the true uncured = 1; cured = 0 status), and
#' covariates. Variables prefixed with \code{"U"} indicates unpenalized
#' covariates and is equal to the value
#' passed to \code{nonp} (default is 2). Variables prefixed with \code{"X"}
#' indicates penalized covariates and is equal to the value passed to
#' \code{j}.
#' }
#' @return \item{testing_y}{the true status for the testing set: uncured = 1;
#' cured = 0}
#' @return \item{parameters}{a list including: the indices of true incidence
#' signals (\code{nonzero_b}), indices of true latency signals
#' (\code{nonzero_beta}), unpenalized incidence parameter values (\code{b_u}),
#' unpenalized latency parameter values (\code{beta_u}), parameter values for
#' the true incidence signals among penalized covariates (\code{b_p_nz}),
#' parameter values for the true latency signals among penalized covariates
#' (\code{beta_p_nz}), parameter value for the incidence intercept
#' (\code{itct})}
#'
#' @export
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstats {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*
#' @import mvnfast
#' @import flexsurv
#' @import stats
#'
#' @examples
#' library(survival)
#' withr::local_seed(1234)
#' # This dataset has 2 penalized features associated with the outcome,
#' # 3 penalized features not associated with the outcome (noise features), and 1
#' # unpenalized noise feature.
#' data <- generate_cure_data(n = 1000, j = 5, n_true = 2, nonp = 1, a = 2)
#' # Extract the training data
#' training <- data$training
#' # Extract the testing data
#' testing <- data$testing
#' # To identify the features truly associated with incidence
#' names(training)[grep("^X", names(training))][data$parameters$nonzero_b]
#' # To identify the features truly associated with latency
#' names(training)[grep("^X", names(training))][data$parameters$nonzero_beta]
#' # Fit the model to the training data
#' fitem <- cureem(Surv(Time, Censor) ~ ., data = training,
#'   x_latency = training)
#' # Examine the estimated coefficients at the (default) minimum AIC
#' coef(fitem)
#' # As the penalty increases, the coefficients for the noise variables shrink
#' # to or remain at zero, while the truly associated features have coefficient
#' # paths that depart from zero. This shows the model's ability to distinguish
#' # signal from noise.
#' plot(fitem, label = TRUE)
generate_cure_data <- function(n = 400, j = 500, nonp = 2, train_prop = 0.75,
                               n_true = 10, a = 1, rho = 0.5, itct_mean = 0.5,
                               cens_ub = 20, alpha = 1, lambda = 2,
                               same_signs = FALSE, model = "weibull") {
  # j: number of penalized covariates
  # nonp: number of non-penalized covariates
  if (n_true > j) {
    stop("Error: \"n_true\" must be less than or equal to j")
  }
  tr_i <- 1:round(train_prop * n) # training index
  te_i <- (round(train_prop * n) + 1):n # testing index  ##corrected Mar30

  ## covariance matrices
  sd <- 0.5
  block_sz <- round(j / n_true)

  corr_x_p <- matrix(0, j, j)
  for (i in 1:min(ceiling(n_true/block_sz+1), n_true)) {
    if (j %% n_true == 0) {
      corr_x_p[(block_sz * (i - 1) + 1):(block_sz * i), (block_sz * (i - 1) + 1):
                 (block_sz * i)] <- rho^abs(outer(1:block_sz, 1:block_sz, "-"))
    } else {
        if (j - block_sz*i > 0) {
          corr_x_p[(block_sz * (i - 1) + 1):(block_sz * i), (block_sz * (i - 1) + 1):
                 (block_sz * i)] <- rho^abs(outer(1:block_sz, 1:block_sz, "-"))
        }
    }
    diag(corr_x_p) <- 1
    }
  sigma_x_p <- sd^2 * corr_x_p
  x_p <- mvnfast::rmvn(n, mu = rep(0, j), sigma = sigma_x_p)
  ## unpenalized
  if (nonp > 0) {
    x_u <- matrix(rnorm(n * nonp, sd = sd), ncol = nonp)
    colnames(x_u) <- paste0("U", 1:dim(x_u)[2])
    w_u <- x_u
    if (!same_signs) { # true signals from two parts are randomly generated
      b_u <- rnorm(nonp, mean = 0.3, sd = 0.1)
      b_u <- b_u * sample(c(1, -1), length(b_u), replace = TRUE)
      beta_u <- rnorm(nonp, mean = 0.3, sd = 0.1)
      beta_u <- beta_u * sample(c(1, -1), length(beta_u), replace = TRUE)
    } else { # true signals from two parts have same signs
      sign_u <- sample(c(1, -1), nonp, replace = TRUE)
      b_u <- abs(rnorm(nonp, mean = 0.3, sd = 0.1)) * sign_u
      beta_u <- abs(rnorm(nonp, mean = 0.3, sd = 0.1)) * sign_u
    }
  }
  w_p <- x_p
  ## penalized
  if (!same_signs) { # true signals from two parts are randomly generated
    nonzero_b <- nonzero_beta <- rep(NA, n_true)
    for (i in 1:n_true) {
      if (i == n_true) {
        if (length(((n_true-1)*block_sz+1):j) == 1) {
          nonzero_b[i] <- j
          nonzero_beta[i] <- j
        } else {
          nonzero_b[i] <- sample(x=((n_true-1)*block_sz+1):j, size=1)
          nonzero_beta[i] <- sample(x=((n_true-1)*block_sz+1):j, size=1)
        }
      } else {
        nonzero_b[i] <- sample((block_sz * (i - 1) + 1):(block_sz * i), 1)
        nonzero_beta[i] <- sample((block_sz * (i - 1) + 1):(block_sz * i), 1)
      }
    }
    b_p <- beta_p <- rep(0, j)
    b_p[nonzero_b] <- a * sample(c(1, -1), n_true, replace = TRUE)
    beta_p[nonzero_beta] <- a * sample(c(1, -1), n_true, replace = TRUE)
  } else { # true signals from two parts have same signs
    nonzero_b <- rep(NA, n_true)
    for (i in 1:n_true) {
      nonzero_b[i] <- sample((block_sz * (i - 1) + 1):
      (block_sz * i), 1)
    }
    b_p <- rep(0, j)
    b_p[nonzero_b] <- a * sample(c(1, -1), n_true, replace = TRUE)
    beta_p <- b_p
    nonzero_beta <- nonzero_b
  }
  ## cure or not
  itct <- rnorm(1, mean = itct_mean, sd = 0.1)
  if (nonp > 0) {
    bx <- itct + x_u %*% b_u + x_p %*% b_p
  } else {
    bx <- itct + x_p %*% b_p
  }
  p <- 1 / (1 + exp(-bx))
  y <- rbinom(n, 1, p)

  ## survival
  if (nonp > 0) {
    beta_w <- w_u %*% beta_u + w_p %*% beta_p
  } else {
    beta_w <- w_p %*% beta_p
  }
  if (model == "weibull") {
    t <- rweibull(n, shape = alpha, scale = 1 / lambda * exp(-beta_w / alpha))
  } else if (model == "GG") {
    t <- flexsurv::rgengamma(n,
      mu = -log(lambda) - beta_w / alpha,
      sigma = 1 / alpha, Q = 2
    )
  } else if (model == "Gompertz") {
    t <- flexsurv::rgompertz(n, shape = 0.2, rate = exp(beta_w))
  } else if (model == "nonparametric") {
    t <- nonparametric_time_generator(n, beta_w, maxT = cens_ub, knots = 8)
  } else if (model == "GG_baseline") {
    t <- rep(NA, n)
    for (i in 1:n) {
      u <- runif(1)
      t[i] <- uniroot(function(a) flexsurv::pgengamma(a, mu = -log(2), sigma = 1, Q = 0.5, lower.tail = FALSE)^exp(beta_w[i]) - u,
        c(0, 20),
        extendInt = "yes"
      )$root
    }
  }
  u <- runif(n, 0, cens_ub)
  delta <- ifelse(t > u | y == 0, 0, 1)
  time <- pmin(t, u)
  time[y == 0] <- u[y == 0]
  ## training and test
  if (nonp > 0) {
    colnames(x_u) <- paste0("U",1:dim(x_u)[2])
    colnames(w_u) <- paste0("U",1:dim(w_u)[2])
    tr_data <- list(
      x_u = x_u[tr_i, , drop=FALSE], x_p = x_p[tr_i, , drop = FALSE],
      w_u = w_u[tr_i, , drop = FALSE],
      w_p = w_p[tr_i, , drop = FALSE], time = time[tr_i], y = y[tr_i],
      delta = delta[tr_i],
      y = y[tr_i]
    )
    te_data <- list(
     x_u = x_u[te_i, , drop = FALSE], x_p = x_p[te_i, , drop = FALSE],
     w_u = w_u[te_i, , drop = FALSE],
     w_p = w_p[te_i, , drop = FALSE], time = time[te_i], y = y[te_i],
     delta = delta[te_i],
     y = y[tr_i]
    )
   training <- data.frame(
      Time = tr_data$time, Censor = tr_data$delta,
      tr_data$x_u, tr_data$x_p
    )
   training_y <- tr_data$y
   testing <- data.frame(
      Time = te_data$time, Censor = te_data$delta,
      te_data$x_u, te_data$x_p
    )
   testing_y <- te_data$y
    parameters <- list(
      nonzero_b = nonzero_b, nonzero_beta = nonzero_beta,
      b_u = b_u, beta_u = beta_u, b_p_nz = b_p[nonzero_b],
      beta_p_nz = beta_p[nonzero_beta], itct = itct
    )
  } else {
    tr_data <- list(
      x_p = x_p[tr_i, , drop = FALSE],
      w_p = w_p[tr_i, , drop = FALSE ], time = time[tr_i], y = y[tr_i],
      delta = delta[tr_i], y = y[tr_i]
    )
    te_data <- list(
      x_p = x_p[te_i, , drop = FALSE],
      w_p = w_p[te_i, , drop = FALSE], time = time[te_i], y = y[te_i],
      delta = delta[te_i], y = y[te_i]
    )
    training <- data.frame(
      Time = tr_data$time, Censor = tr_data$delta,
      tr_data$x_p
    )
    training_y <- tr_data$y
    testing <- data.frame(
      Time = te_data$time, Censor = te_data$delta,
      te_data$x_p
    )
    testing_y <- te_data$y
    parameters <- list(
      nonzero_b = nonzero_b, nonzero_beta = nonzero_beta,
      b_p_nz = b_p[nonzero_b],
      beta_p_nz = beta_p[nonzero_beta], itct = itct
    )
  }
  return(list(training = training, training_y = training_y, testing = testing,
              testing_y = testing_y, parameters = parameters))
}
