#' Non-parametric pest for a non-zero cured fraction
#'
#' @description
#' Tests the null hypothesis that the proportion of observations susceptible to
#' the event = 1 against the alternative that the proportion of observations
#' susceptible to the event is < 1. If the null hypothesis is rejected, there
#' is a significant cured fraction.
#'
#' @param object  a \code{survfit} object.
#' @param reps number of simulations on which to base the p-value
#' (default = 1000).
#' @param seed optional random seed.
#' @param plot logical. If TRUE a histogram of the estimated susceptible
#' proportions over all simulations is produced.
#' @param b optional. If specified the maximum observed time for the uniform
#' distribution for generating the censoring times. If not specified, an
#' exponential model is used for generating the censoring times (default).
#'
#' @return \item{proportion_susceptible}{estimated proportion of susceptibles}
#' @return \item{proportion_cured}{estimated proportion of those cured}
#' @return \item{p_value}{p-value testing the null hypothesis that the
#' proportion of susceptibles = 1 (cured fraction = 0) against the alternative
#' that the proportion of susceptibles < 1 (non-zero cured fraction)}
#' @return \item{time_95_percent_of_events}{estimated time at which 95% of
#' events should have occurred}
#' @export
#'
#' @references Maller, R. A. and Zhou, X. (1996) \emph{Survival Analysis with
#' Long-Term Survivors}. John Wiley & Sons.
#'
#' @srrstats {G1.0} *Statistical Software should list at least one primary reference from published academic literature.*
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
#' @srrstats {G5.2} *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstats {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*
#' @srrstats {RE4.5} *Numbers of observations submitted to model (via `nobs()`)*
#' @seealso \code{\link{survfit}}, \code{\link{cure_estimate}},
#' \code{\link{sufficient_fu_test}}
#'
#' @import survival
#' @import stats
#' @import flexsurv
#' @import flexsurvcure
#' @keywords htest
#' @examples
#' library(survival)
#' withr::local_seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' km_fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
#' nonzerocure_test(km_fit)
nonzerocure_test <- function(object, reps = 1000, seed = NULL, plot = FALSE,
                             b = NULL) {
  if (!(c("survfit") %in% class(object))) {
    stop("Error: object must be a survfit object")
  }
  # Surv(Time, Censor)
  Y <- Surv(object$time, object$n.event)
  n <- length(Y)
  if (is.null(object$strata)) {
    exp_model <- survival::survreg(Y ~ 1,
      data = Y,
      dist = "exponential"
    )
    mcm <- flexsurvcure::flexsurvcure(Y ~ 1,
      data = Y,
      dist = "exp", mixture = TRUE
    )
    susceptible <- 1 - cure_estimate(object)
    mean <- unname(exp(exp_model$coefficients["(Intercept)"]))
    rate <- unname(mcm$coefficients["rate"])
    cured <- cure_estimate(object)
    if (cured > 0 && cured < 1) {
      time_95 <- 1 / exp(rate) - log((0.05 * cured) / (0.95 * (1 - cured)))
    } else {
      time_95 <- 1 / exp(rate)
    }
    p <- sim_cure(n,
      mu = mean, censor_mu = time_95, reps = reps, seed = seed,
      b = b
    )
    p_value <- ifelse(sum(p < susceptible) / reps == 0, paste("<", 1 / reps),
      sum(p < susceptible) / reps
    )
    if (plot) {
      hist(p,
        xlab = "Proportion Susceptible",
        main = "Non-parametric Reference Distribution"
      )
    }
  } else {
    strata <- gsub(".*=", "", names(object$strata))
    mean <- rep(NA, length(strata))
    rate <- rep(NA, length(strata))
    for (i in seq_along(strata)) {
      group <- extract_rhs_values(object)
      frame <- data.frame(surv = Y, group = group)
      exp_model <- survival::survreg(surv ~ 1,
        dist = "exponential", data = frame,
        subset = (group == strata[i]))
      mean[i] <- exp(exp_model$coefficients["(Intercept)"])
      mcm <- flexsurvcure::flexsurvcure(surv ~ 1,
        data = frame, subset = group == strata[i],
        dist = "exp", mixture = TRUE
      )
      rate[i] <- mcm$coefficients["rate"]
      # }
    }
    mean <- unname(mean)
    rate <- unname(rate)
    cured <- cure_estimate(object)[, "Cured"]
    susceptible <- cure_estimate(object)
    n <- unname(object$strata)
    time_95 <- numeric()
    for (i in seq_along(strata)) {
      if (cured[i] > 0 && cured[i] < 1) {
        time_95[i] <- 1 / exp(rate[i]) - log((0.05 * cured[i]) /
          (0.95 * (1 - cured[i])))
      } else {
        time_95[i] <- 1 / exp(rate[i])
      }
    }
    p <- list()
    for (i in seq_along(object$strata)) {
      p[[i]] <- sim_cure(n[i],
        mu = mean[i], censor_mu = time_95[i],
        reps = reps, seed = seed, b = b
      )
    }
    p_value <- numeric()
    for (i in seq_along(object$strata)) {
      p_value[i] <- ifelse(sum(p[[i]] <
        susceptible[i, "Susceptible"]) / reps == 0,
      paste("<", 1 / reps),
      sum(p[[i]] < susceptible[i, "Susceptible"]) / reps
      )
    }
    names(p_value) <- names(object$strata)
    if (plot) {
      for (i in seq_along(p)) {
        hist(p[[i]],
          xlab = "Proportion Susceptible",
          main = paste(names(object$strata[i]),
            "Non-parametric Reference Distribution",
            sep = " "
          )
        )
      }
    }
    susceptible <- susceptible[, "Susceptible"]
    names(time_95) <- names(susceptible) <- names(cured) <- names(object$strata)
  }
  list(
    proportion_susceptible = susceptible, proportion_cured = cured,
    p_value = p_value, time_95_percent_of_events = time_95
  )
}
