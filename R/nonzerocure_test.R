#' Non-parametric pest for a non-zero cured fraction
#'
#' @description
#' Tests the null hypothesis that the proportion of observations susceptible to the event = 1 against the alternative that the proportion of observations susceptible to the event is < 1. If the null hypothesis is rejected, there is a significant cured fraction.
#'
#' @param object  a \code{survfit} object.
#' @param Reps number of simulations on which to base the p-value (default = 1000).
#' @param seed optional random seed.
#' @param plot logical. If TRUE a histogram of the estimated susceptible proportions over all simulations is produced.
#' @param B optional. If specified the maximum observed time for the uniform distribution for generating the censoring times. If not specified, an exponential model is used for generating the censoring times (default).
#'
#' @return \item{proportion_susceptible}{estimated proportion of susceptibles}
#' @return \item{proportion_cured}{estimated proportion of those cured}
#' @return \item{p.value}{p-value testing the null hypothesis that the proportion of susceptibles = 1 (cured fraction = 0) against the alternative that the proportion of susceptibles < 1 (non-zero cured fraction)}
#' @return \item{time_95_percent_of_events}{estimated time at which 95% of events should have occurred}
#' @export
#'
#' @references Maller, R. A. and Zhou, X. (1996) \emph{Survival Analysis with Long-Term Survivors}. John Wiley & Sons.
#'
#' @seealso \code{\link{survfit}}, \code{\link{cure_estimate}}, \code{\link{sufficient_fu_test}}
#'
#' @import survival
#' @import stats
#' @import flexsurv
#' @import flexsurvcure
#' @keywords htest
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(N = 100, J = 10, nTrue = 10, A = 1.8)
#' training <- temp$Training
#' km.fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
#' nonzerocure_test(km.fit)
nonzerocure_test <- function(object, Reps = 1000, seed = NULL, plot =FALSE, B = NULL) {
  if (!(c("survfit") %in% class(object))) stop("object must be a survfit object")
  if (is.null(object$strata)) {
    exp.model <- survival::survreg(as.formula(object$call$formula), data=eval(parse(text=object$call$data)), dist="exponential")
    mcm <- flexsurvcure::flexsurvcure(as.formula(object$call$formula), data=eval(parse(text=object$call$data)), dist="exp", mixture=T)
    susceptible<-1-cure_estimate(object)
    mean <- unname(exp(exp.model$coefficients["(Intercept)"]))
    rate <- unname(mcm$coefficients["rate"])
    cured <- cure_estimate(object)
    if (class(eval(parse(text=object$call$data))) %in% c("ExpressionSet", "SummarizedExperiment")) {
      n <- dim(eval(parse(text=object$call$data)))[2]
    } else {
      n <- dim(eval(parse(text=object$call$data)))[1]
    }
    if (cured>0 & cured <1) {
      time.95 <- 1/exp(rate)-log((0.05*cured)/(0.95*(1-cured)))
    } else {
      time.95 <- 1/exp(rate)
    }
    p <- sim_cure(n, mu=mean, censor.mu=time.95, Reps = Reps, seed=seed, B=B)
    p.value<-ifelse(sum(p < susceptible)/Reps==0,paste("<",1/Reps),sum(p<susceptible)/Reps)
    if (plot) {
      hist(p, xlab="Proportion Susceptible", main="Non-parametric Reference Distribution")
    }
  } else {
    strata <- gsub(".*=","",names(object$strata))
    mean <- numeric()
    rate <- numeric()
    for (i in 1:length(strata)) {
      ff <- as.formula(object$call$formula)
      tt <- terms(ff)
      tt[-1]
      rhs.variable <- all.vars(as.formula(object$call$formula))[-(1:2)]
      if (class(eval(parse(text=object$call$data))) %in% c("ExpressionSet", "SummarizedExperiment")) {
        exp.model <- survival::survreg(tt[-1], data=(eval(parse(text=object$call$data))), dist="exponential", subset=((eval(parse(text=object$call$data)))[,rhs.variable]==strata[i]))
        mean <- c(mean, exp(exp.model$coefficients["(Intercept)"]))
        mcm <- flexsurvcure::flexsurvcure(formula(tt[-1]), data=(eval(parse(text=object$call$data))), subset=((eval(parse(text=object$call$data)))[,rhs.variable]==strata[i]), dist="exp", mixture=T)
        rate <- c(rate, mcm$coefficients["rate"])
      }
    }
    mean <- unname(mean)
    rate <- unname(rate)
    cured <- cure_estimate(object)[,"Cured"]
    susceptible<-cure_estimate(object)
    n <- unname(object$strata)
    time.95 <- numeric()
    for (i in 1:length(strata)) {
      if (cured[i]>0 & cured[i] <1) {
        time.95[i] <- 1/exp(rate[i])-log((0.05*cured[i])/(0.95*(1-cured[i])))
      } else {
        time.95[i] <- 1/exp(rate[i])
      }
    }
    p <- list()
    for (i in 1:length(object$strata)) {
      p[[i]] <- sim_cure(n[i], mu=mean[i], censor.mu=time.95[i], Reps = Reps, seed=seed, B=B)
    }
    p.value<-numeric()
    for (i in 1:length(object$strata)) {
      p.value[i]<-ifelse(sum(p[[i]] < susceptible[i,"Susceptible"])/Reps==0,paste("<",1/Reps),sum(p[[i]]<susceptible[i,"Susceptible"])/Reps)
    }
    names(p.value)<-names(object$strata)
    if (plot) {
      for (i in 1:length(p)) {
        hist(p[[i]], xlab="Proportion Susceptible", main=paste(names(object$strata[i]),"Non-parametric Reference Distribution", sep=" "))
      }
    }
    susceptible<-susceptible[,"Susceptible"]
    names(time.95) <- names(susceptible) <- names(cured) <- names(object$strata)
  }
  list(proportion_susceptible=susceptible, proportion_cured = cured, p.value=p.value, time_95_percent_of_events = time.95)
}
