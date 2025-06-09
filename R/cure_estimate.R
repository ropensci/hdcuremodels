#' Estimate cured fraction
#'
#' @description
#' Estimates the cured fraction using a Kaplan-Meier fitted object.
#'
#' @param object a \code{survfit} object.
#'
#' @return estimated proportion of cured observations
#' @export
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {G1.3} *All statistical terminology should be clarified and unambiguously defined.*
#' @srrstats {G2.0} *Implement assertions on lengths of inputs, particularly through asserting that inputs expected to be single- or multi-valued are indeed so.*
#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @seealso \code{\link[survival]{survfit}}, \code{\link{sufficient_fu_test}},
#' \code{\link{nonzerocure_test}}
#' @import survival
#' @import stats
#' @keywords univar
#'
#' @examples
#' library(survival)
#' withr::local_seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' km_fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
#' cure_estimate(km_fit)
cure_estimate <- function(object) {
  if (!(c("survfit") %in% class(object))) {
    stop("Error: object must be a survfit object")
  }
  summary_kme <- summary(object)
  is_strata <- grep("strata", names(summary_kme))
  if (length(is_strata) != 0) {
    min_surv <- aggregate(summary_kme$surv,
      by = list(summary_kme$strata),
      min
    )
    colnames(min_surv)[2] <- "Cured"
    min_surv$"Susceptible" <- 1 - min_surv$Cured
  } else {
    min_surv <- min(summary_kme[["surv"]])
  }
  min_surv
}
