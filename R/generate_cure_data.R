#' Simulate data under a mixture cure model
#'
#' @param N an integer denoting the total sample size.
#' @param J an integer denoting the number of penalized predictors which is the same for both the incidence and latency portions of the model.
#' @param nonp an integer less than J denoting the number of unpenalized predictors (which is the same for both the incidence and latency portions of the model.
#' @param train.prop a numeric value in 0, 1 representing the fraction of N to be used in forming the Training dataset.
#' @param nTrue an integer denoting the number of variables truly associated with the outcome (i.e., the number of covariates with nonzero parameter values) among the penalized predictors.
#' @param A a numeric value denoting the effect size which is the same for both the incidence and latency portions of the model.
#' @param rho a numeric value in 0, 1 representing the correlation between adjacent covariates in the same block. See details below.
#' @param itct_mean a numeric value representing the expectation of the incidence intercept which controls the cure rate.
#' @param cens_ub a numeric value representing the upper bound on the censoring time distribition which follows a uniform distribution on 0, \code{cens_ub}.
#' @param alpha a numeric value representing the shape parameter in the Weibull density.
#' @param lambda a numeric value representing the rate parameter in the Weibull density.
#' @param same_signs logical, if TRUE the incidence and latency coefficients have the same signs.
#' @param model type of regression model to use for the latency portion of mixture cure model. Can be "weibull", "GG", "Gompertz", "nonparametric", or "GG_baseline".
#'
#' @return \item{Training}{Training data.frame which includes Time, Censor, and covariates.}
#' @return \item{Testing}{Testing data.frame which includes Time, Censor, and covariates.}
#' @return \item{parameters}{A list including: the indices of true incidence signals (\code{nonzero_b}), indices of true latency signals (\code{nonzero_beta}), unpenalized incidence parameter values (\code{b_u}), unpenalized latency parameter values (\code{beta_u}), parameter values for the true incidence signals among penalized covariates (\code{b_p_nz}), parameter values for the true latency signals among penalized covariates (\code{beta_p_nz}), parameter value for the incidence intercept (\code{itct})}
#'
#' @export
#'
#' @import mvnfast
#' @import flexsurv
#' @import stats
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' data <- generate_cure_data(N = 200, J = 50, nTrue = 10, A = 1.8, rho = 0.2)
#' training <- data$Training
#' testing <- data$Testing
#' fit <- cureem(Surv(Time, Censor) ~ ., data = training,
#'               x.latency = training, model = "cox", penalty = "lasso",
#'               lambda.inc = 0.05, lambda.lat = 0.05,
#'               gamma.inc = 6, gamma.lat = 10)
generate_cure_data <- function(N=400, J=500, nonp=2, train.prop = 3/4, nTrue=10, A=1, rho=0.5, itct_mean=0.5, cens_ub = 20,
         alpha=1, lambda=2, same_signs=FALSE, model="weibull"){
  # J: number of penalized covariates
  # nonp: number of non-penalized covariates
  tr_i = 1:round(train.prop*N) # training index
  te_i = (round(train.prop*N)+1):N # testing index  ##corrected Mar30

  ## covariance matrices
  sd = 0.5
  block_sz = round(J/nTrue)
  corr_X_p = matrix(0, J, J)
  for (i in 1:nTrue)
    corr_X_p[(block_sz*(i-1)+1):(block_sz*i),(block_sz*(i-1)+1):(block_sz*i)] =
    rho^abs(outer(1:block_sz, 1:block_sz, "-"))
  Sigma_X_p = sd^2*corr_X_p
  X_p = mvnfast::rmvn(N, mu = rep(0,J), sigma = Sigma_X_p)
  X_u = matrix(rnorm(N*nonp, sd=sd), ncol = nonp)
  W_p = X_p
  W_u = X_u

  if(!same_signs){ # true signals from two parts are randomly generated
    ## unpenalized
    b_u = rnorm(nonp, mean = 0.3, sd = 0.1)
    b_u = b_u * sample(c(1,-1), length(b_u), replace = T)
    beta_u = rnorm(nonp, mean = 0.3, sd = 0.1)
    beta_u = beta_u * sample(c(1,-1), length(beta_u), replace = T)

    ## penalized
    nonzero_b <- nonzero_beta <- rep(NA, nTrue)
    for(i in 1:nTrue){
      nonzero_b[i] = sample((block_sz*(i-1)+1):(block_sz*i), 1)
      nonzero_beta[i] = sample((block_sz*(i-1)+1):(block_sz*i), 1)
    }
    b_p <- beta_p <- rep(0, J)
    b_p[nonzero_b] = A * sample(c(1,-1), nTrue, replace = T)
    beta_p[nonzero_beta] = A * sample(c(1,-1), nTrue, replace = T)
  } else { # true signals from two parts have same signs
    ## unpenalized
    sign_u = sample(c(1,-1), nonp, replace = T)
    b_u = abs(rnorm(nonp, mean = 0.3, sd = 0.1)) * sign_u
    beta_u = abs(rnorm(nonp, mean = 0.3, sd = 0.1)) * sign_u
    ## penalized
    nonzero_b <- rep(NA, nTrue)
    for(i in 1:nTrue) nonzero_b[i] = sample((block_sz*(i-1)+1):(block_sz*i), 1)
    b_p  <- rep(0, J)
    b_p[nonzero_b] = A * sample(c(1,-1), nTrue, replace = T)
    beta_p = b_p
    nonzero_beta = nonzero_b
  }

  ## cure or not
  itct = rnorm(1, mean = itct_mean, sd = 0.1)
  bx = itct + X_u %*% b_u + X_p %*% b_p
  p = 1/(1+exp(-bx))
  Y = rbinom(N, 1, p)

  ## survival
  beta_w = W_u%*%beta_u + W_p%*%beta_p
  if (model=="weibull") {
    t = rweibull(N, shape = alpha, scale = 1/lambda * exp(- beta_w/alpha))
  } else if (model=="GG") {
    t = flexsurv::rgengamma(N, mu=-log(lambda)-beta_w/alpha, sigma=1/alpha,
                            Q=2)
  } else if (model=="Gompertz") {
    t = flexsurv::rgompertz(N, shape = 0.2, rate = exp(beta_w))
  } else if (model=="nonparametric") {
    t = nonparametric_time_generator(N, beta_w, maxT = cens_ub, knots = 8)
  } else if (model=="GG_baseline") {
    t = rep(NA, N)
    for(i in 1:N){
      u = runif(1)
      t[i] = uniroot(function(a)
        flexsurv::pgengamma(a, mu=-log(2), sigma=1, Q=0.5, lower.tail = FALSE)^exp(beta_w[i])-u,
        c(0,20), extendInt = "yes")$root
    }
  }
  u = runif(N, 0, cens_ub)
  delta = ifelse(t>u | Y==0, 0, 1)
  time = pmin(t, u)
  time[Y==0] = u[Y==0]

  ## training and test
  Tr_data = list(X_u=X_u[tr_i,], X_p=X_p[tr_i,], W_u=W_u[tr_i,], W_p=W_p[tr_i,], time=time[tr_i], Y=Y[tr_i],
                 delta=delta[tr_i])
  Te_data = list(X_u=X_u[te_i,], X_p=X_p[te_i,], W_u=W_u[te_i,], W_p=W_p[te_i,], time=time[te_i], Y=Y[te_i], delta=delta[te_i])
  Training <- data.frame(Time = Tr_data$time, Censor = Tr_data$delta, Tr_data$X_u, Tr_data$X_p)
  Testing <- data.frame(Time = Te_data$time, Censor = Te_data$delta, Te_data$X_u, Te_data$X_p)
  parameters<- list(nonzero_b=nonzero_b, nonzero_beta=nonzero_beta,
                    b_u = b_u, beta_u = beta_u, b_p_nz = b_p[nonzero_b], beta_p_nz = beta_p[nonzero_beta],
                    itct=itct)
  return(list(Training = Training, Testing = Testing, parameters = parameters ))
}
