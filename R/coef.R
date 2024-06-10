#' Extract model coefficients from a fitted mixture cure object
#'
#' @description
#' \code{coef.mixturecure} is a generic function which extracts the model coefficients from a fitted mixture cure model object fit using \code{curegmifs}, \code{cureem}, \code{cv_curegmifs}, or \code{cv_cureem}.
#'
#' @param object a \code{mixturecure} object resulting from \code{curegmifs}, \code{cureem}, \code{cv_curegmifs}, or \code{cv_cureem}.
#' @param model.select for models fit using \code{curegmifs} or \code{cureem} any step along the solution path can be selected. The default is \code{model.select = "AIC"} which calculates the predicted values using the coefficients from the model having the lowest AIC. Other options are \code{model.select = "mAIC"} for the modified AIC,  \code{model.select = "cAIC"} for the corrected AIC, \code{model.select = "BIC"}, \code{model.select = "mBIC"} for the modified BIC, \code{model.select = "EBIC"} for the extended BIC, \code{model.select = "logLik"} for the step that maximizes the log-likelihood, or any numeric value from the solution path. This option has no effect for objects fit using \code{cv_curegmifs} or \code{cv_cureem}.
#' @param ... other arguments.
#'
#' @return a list of estimated parameters extracted from the model object using the model selection criterion
#' @seealso \code{\link{curegmifs}}, \code{\link{cureem}}, \code{\link{summary.mixturecure}}, \code{\link{plot.mixturecure}}, \code{\link{predict.mixturecure}}
#' @export
#' @keywords methods
#' @method coef mixturecure
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(N = 100, J = 10, nTrue = 10, A = 1.8)
#' training <- temp$Training
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'                       data = training, x.latency = training,
#'                       model = "weibull", thresh = 1e-4, maxit = 2000, epsilon = 0.01,
#'                       verbose = FALSE)
#' coef(fit)
coef.mixturecure<-
  function (object, model.select = "AIC", ...)
  {
    if (!(class(object) %in% "mixturecure")) stop("class of object must be mixturecure")
    if (!object$cv) {
      if (is.character(model.select)) {
        model.select <- c("AIC", "BIC", "logLik", "cAIC", "mAIC", "mBIC", "EBIC")[pmatch(model.select,
                                                                                         c("AIC", "BIC", "logLik", "cAIC", "mAIC", "mBIC", "EBIC"))]
        if (any(!model.select%in%c("AIC", "BIC", "logLik", "cAIC", "mAIC", "mBIC", "EBIC")))
          stop("model.select must be either 'AIC', 'BIC', 'logLik', 'cAIC', 'mAIC', 'mBIC', or 'EBIC' ")
        if (!is.null(object$x.incidence)) {
          vars.inc <- apply(object$b_path, 1, function(x) sum(x !=
                                                              0))
        }
        else {
          vars.inc <- 0
        }
        if (!is.null(object$x.latency)) {
          vars.lat <- apply(object$beta_path, 1, function(x) sum(x !=
                                                                 0))
        }
        else {
          vars.lat <- 0
        }
      }
      else {
        if (!is.null(object$x.latency)) {
          vars.lat <- sum(object$beta != 0)
        }
        else {
          vars.lat <- 0
        }
        if (!is.null(object$x.incidence)) {
          vars.inc <- sum(object$b != 0)
        }
        else {
          vars.inc <- 0
        }
      }
      if (object$model == "weibull") {
        df <- vars.inc + vars.lat + 3
      }
      else if (object$model == "exponential") {
        df <- vars.inc + vars.lat + 2
      }
      else if (object$model == "cox") {
        df <- vars.inc + vars.lat + 1
      }
      if (object$method == "EM") {
        logLik <- object$logLik.inc + object$logLik.lat
      }
      else {
        logLik <- object$logLik
      }
      p <- dim(object$x.incidence)[2] + dim(object$x.latency)[2]
      AIC <- 2 * df - 2 * logLik
      #cAIC <- AIC+(2*df^2+6*df+4)/(length(object$y)-df-2)
      cAIC<-AIC+(2*df*(df+1))/(length(object$y)-df-1) # https://www.mathworks.com/help/econ/information-criteria.html
      mAIC <- (2+2*log(p/.5)) * df - 2 * logLik
      BIC <- df * (log(length(object$y))) -  2 * logLik
      mBIC <- df * (log(length(object$y)) + 2*log(p/4)) - 2 * logLik
      EBIC <- log(length(object$y)) * df + 2*(1-.5)*log(choose(p, df)) - 2 * logLik
    if (object$model != "cox") {
      if (object$mode == "weibull") {
        if (!exists("alpha_path", object))
        object$alpha_path <- object$alpha
      }
      if (!exists("rate_path", object))
        object$rate_path <- object$rate
    }
      if (!object$cv) {
        if (model.select == "AIC") {
          model.select = which.min(AIC)
        }
        else if (model.select == "BIC") {
          model.select = which.min(BIC)
        }
        else if (model.select == "mAIC") {
          model.select = which.min(mAIC)
        }
        else if (model.select == "mBIC") {
          model.select = which.min(mBIC)
        }
        else if (model.select == "EBIC") {
          model.select = which.min(EBIC)
        }
        else if (model.select == "cAIC") {
          model.select = which.min(cAIC)
        }
        else if (model.select == "logLik") {
          model.select = which.max(logLik)
        }
      }
      if (is.null(object$x.latency) & !is.null(object$x.incidence)) {
        if (object$model == "exponential") {
          coef <- list(rate = object$rate_path[model.select],
                       b0 = object$b0_path[model.select], beta_inc = object$b_path[model.select,
                       ])
        }
        if (object$model == "weibull") {
          coef <- list(rate = object$rate_path[model.select],
                       shape = object$alpha_path[model.select], b0 = object$b0_path[model.select],
                       beta_inc = object$b_path[model.select, ])
        }
        if (object$model == "cox") {
          coef <- list(b0 = object$b0_path[model.select],
                       beta_inc = object$b_path[model.select, ])
        }
      }
      else if (!is.null(object$x.latency) & is.null(object$x.incidence)) {
        if (object$model == "exponential") {
          coef <- list(rate = object$rate_path[model.select],
                       beta_lat = object$beta_path[model.select, ])
        }
        if (object$model == "weibull") {
          coef <- list(rate = object$rate_path[model.select],
                       shape = object$alpha_path[model.select], beta_lat = object$beta_path[model.select,
                       ])
        }
        if (object$model == "cox") {
          coef <- list(beta_lat = object$beta_path[model.select,
          ])
        }
      }
      else if (!is.null(object$x.latency) & !is.null(object$x.incidence)) {
        if (object$model == "exponential") {
          coef <- list(rate = object$rate_path[model.select],
                       b0 = object$b0_path[model.select], beta_inc = object$b_path[model.select,
                       ], beta_lat = object$beta_path[model.select,
                       ])
        }
        if (object$model == "weibull") {
          coef <- list(rate = object$rate_path[model.select],
                       shape = object$alpha_path[model.select], b0 = object$b0_path[model.select],
                       beta_inc = object$b_path[model.select, ], beta_lat = object$beta_path[model.select,
                       ])
        }
        if (object$model == "cox") {
          coef <- list(b0 = object$b0_path[model.select],
                       beta_inc = object$b_path[model.select, ], beta_lat = object$beta_path[model.select,
                       ])
        }
      }
    }
    else {
      if (object$model == "weibull") {
        if (!is.null(object$x.incidence) & !is.null(object$x.latency))
          coef <- list(rate = object$rate, shape = object$alpha,
                       b0 = object$b0, beta_inc = object$b, beta_lat = object$beta)
        if (is.null(object$x.incidence) & !is.null(object$x.latency))
          coef <- list(rate = object$rate, shape = object$alpha,
                       beta_lat = object$beta)
        if (!is.null(object$x.incidence) & is.null(object$x.latency))
          coef <- list(rate = object$rate, shape = object$alpha,
                       b0 = object$b0, beta_inc = object$b)
      }
      else if (object$model == "exponential") {
        if (!is.null(object$x.incidence) & !is.null(object$x.latency))
          coef <- list(rate = object$rate, b0 = object$b0,
                       beta_inc = object$b, beta_lat = object$beta)
        if (is.null(object$x.incidence) & !is.null(object$x.latency))
          coef <- list(rate = object$rate, beta_lat = object$beta)
        if (!is.null(object$x.incidence) & is.null(object$x.latency))
          coef <- list(rate = object$rate, b0 = object$b0,
                       beta_inc = object$b)
      }
      else if (object$model == "cox") {
        if (!is.null(object$x.incidence) & !is.null(object$x.latency)) {
          coef <- list(b0 = object$b0, beta_inc = object$b, beta_lat = object$beta)
        }
        if (is.null(object$x.incidence) & !is.null(object$x.latency)) {
          coef <- list(beta_lat = object$beta)
        }
        if (!is.null(object$x.incidence) & is.null(object$x.latency)) {
          coef <- list(b0 = object$b0, beta_inc = object$b)
        }
      }
    }
    coef
  }
