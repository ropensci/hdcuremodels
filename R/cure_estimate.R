#' Estimate cured fraction
#'
#' @description
#' Estimates the cured fraction using a Kaplan-Meier fitted object.
#'
#' @param object a \code{survfit} object.
#'
#' @return estimated proportion of cured observations
#' @export
#'
#' @seealso \code{\link{survfit}}, \code{\link{sufficient_fu_test}},
#' \code{\link{nonzerocure_test}}
#' @import survival
#' @import stats
#' @keywords univar
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' km_fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
#' cure_estimate(km_fit)
cure_estimate <- function(object) {
  if (!(c("survfit") %in% class(object))) {
    stop("object must be a survfit object")
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
