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
#' @seealso \code{\link{survfit}}, \code{\link{sufficient_fu_test}}, \code{\link{nonzerocure_test}}
#' @import survival
#' @import stats
#' @keywords univar
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(N = 100, J = 10, nTrue = 10, A = 1.8)
#' training <- temp$Training
#' km.fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
#' cure_estimate(km.fit)
cure_estimate <-
  function(object) {
    if (!(c("survfit") %in% class(object))) stop("object must be a survfit object")
    summary.kme<-summary(object)
    is.strata<-grep("strata",names(summary.kme))
    if (length(is.strata)!=0) {
      min.surv<-aggregate(summary.kme$surv, by=list(summary.kme$strata), min)
      colnames(min.surv)[2]<-"Cured"
      min.surv$"Susceptible"<-1-min.surv$Cured
    } else {
      min.surv<-min(summary.kme[["surv"]])
    }
    min.surv
  }
