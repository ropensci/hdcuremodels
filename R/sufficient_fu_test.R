#' Test for sufficient follow-up
#'
#' @description
#' Tests for sufficient follow-up using a Kaplan-Meier fitted object.
#'
#' @param object a \code{survfit} object.
#'
#' @return \item{p.value}{p-value from testing the null hypothesis that there was not sufficient follow-up against the alternative that there was sufficient follow-up}
#' @return \item{Nn}{total number of events that occurred at time > pmax(0, 2*(last observed event time)-(last observed time)) and < the last observed event time}
#' @return \item{N}{number of observations in the dataset}
#' @export
#'
#' @references Maller, R. A. and Zhou, X. (1996) \emph{Survival Analysis with Long-Term Survivors}. John Wiley & Sons.
#'
#' @seealso \code{\link{survfit}}, \code{\link{cure_estimate}}, \code{\link{nonzerocure_test}}
#' @import survival
#' @import stats
#' @keywords htest
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(N = 100, J = 10, nTrue = 10, A = 1.8)
#' training <- temp$Training
#' km.fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
#' sufficient_fu_test(km.fit)
sufficient_fu_test<-function(object) {
  if (!(c("survfit") %in% class(object))) stop("object must be a survfit object")
  if (length(object$strata)>1) {
    summary.kme<-summary(object)
    object$strata.group<-rep(names(object$strata),object$strata)
    tn<-aggregate(object$time, by=list(object$strata.group), max)
    colnames(tn)[2]<-"tn"
    tnstar<-aggregate(object$time[object$n.event!=0],by=list(object$strata.group[object$n.event!=0]),max)
    colnames(tnstar)[2]<-"tnstar"
    max.surv<-merge(tn, tnstar, by="Group.1")
    max.surv$delta<-max.surv$tn-max.surv$tnstar
    max.surv$Lower<-2*max.surv$tnstar-max.surv$tn
    Nn<-qn<-alpha_n<-numeric()
    for (i in 1:dim(max.surv)[1]) {
      events<-object$n.event[object$strata.group==names(object$strata)[i]]
      times<-object$time[object$strata.group==names(object$strata)[i]]
      Nn[i]<-sum(events[times>max.surv$Lower[i] & times<=max.surv$tnstar[i]])
      qn[i]<-Nn[i]/object$n[i]
      alpha_n[i]<-(1-Nn[i]/object$n[i])^object$n[i]
    }
    output<-data.frame(alpha_n=alpha_n, Nn=Nn, qn=qn, n=object$n)
    rownames(output)<-names(object$strata)
  } else {
    tn<-max(object$time)
    tnstar<-max(object$time[object$n.event!=0])
    delta<-tn-tnstar
    Lower<-2*tnstar-tn
    Nn<-sum(object$n.event[object$time>Lower & object$time<=tnstar])
    alpha_n <- (1-Nn/object$n)^object$n
    qn<-Nn/object$n
    output<-data.frame(p.value=alpha_n, Nn=Nn, N=object$n)
  }
  output
}
