#' Test for sufficient follow-up
#'
#' @description
#' Tests for sufficient follow-up using a Kaplan-Meier fitted object.
#'
#' @param object a \code{survfit} object.
#'
#' @return \item{p_value}{p-value from testing the null hypothesis that there
#' was not sufficient follow-up against the alternative that there was
#' sufficient follow-up}
#' @return \item{n_n}{total number of events that occurred at time >
#' pmax(0, 2*(last observed event time)-(last observed time)) and < the last
#' observed event time}
#' @return \item{N}{number of observations in the dataset}
#' @export
#'
#' @references Maller, R. A. and Zhou, X. (1996) \emph{Survival Analysis with
#' Long-Term Survivors}. John Wiley & Sons.
#'
#' @seealso \code{\link{survfit}}, \code{\link{cure_estimate}},
#' \code{\link{nonzerocure_test}}
#' @import survival
#' @import stats
#' @keywords htest
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' km_fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
#' sufficient_fu_test(km_fit)
sufficient_fu_test <- function(object) {
  if (!(c("survfit") %in% class(object)))
    stop("object must be a survfit object")
  if (length(object$strata) > 1) {
    object$strata.group <- rep(names(object$strata), object$strata)
    tn <- aggregate(object$time, by = list(object$strata.group), max)
    colnames(tn)[2] <- "tn"
    tnstar <- aggregate(object$time[object$n.event != 0],
                        by = list(object$strata.group[object$n.event != 0]),
                        max)
    colnames(tnstar)[2] <- "tnstar"
    max_surv <- merge(tn, tnstar, by = "Group.1")
    max_surv$delta <- max_surv$tn - max_surv$tnstar
    max_surv$lower <- 2 * max_surv$tnstar - max_surv$tn
    n_n <- q_n <- alpha_n <- numeric()
    for (i in seq_len(dim(max_surv)[1])) {
      events <- object$n.event[object$strata.group == names(object$strata)[i]]
      times <- object$time[object$strata.group == names(object$strata)[i]]
      n_n[i] <- sum(events[times > max_surv$lower[i] & times <=
                             max_surv$tnstar[i]])
      q_n[i] <- n_n[i] / object$n[i]
      alpha_n[i] <- (1 - n_n[i] / object$n[i]) ^ object$n[i]
    }
    output <- data.frame(alpha_n = alpha_n, n_n = n_n, q_n = q_n, n = object$n)
    rownames(output) <- names(object$strata)
  } else {
    tn <- max(object$time)
    tnstar <- max(object$time[object$n.event != 0])
    lower <- 2 * tnstar - tn
    n_n <- sum(object$n.event[object$time > lower & object$time <= tnstar])
    alpha_n <- (1 - n_n / object$n) ^ object$n
    q_n <- n_n / object$n
    output <- data.frame(p_value = alpha_n, n_n = n_n, N = object$n)
  }
  output
}
