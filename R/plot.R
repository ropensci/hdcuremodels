#' Plot fitted mixture cure model
#'
#' @description
#' This function plots either the coefficient path, the AIC, the cAIC, the BIC,
#' or the log-likelihood for a fitted \code{curegmifs} or \code{cureem} object.
#' This function produces a lollipop plot of the coefficient estimates for a
#' fitted \code{cv_curegmifs} or \code{cv_cureem} object.
#'
#' @param x a \code{mixturecure} object resulting from \code{curegmifs} or
#' \code{cureem}, \code{cv_curegmifs} or \code{cv_cureem}.
#' @param type a case-sensitive parameter with default \code{"trace"} which plots the coefficient path for
#' the fitted object. Also available are \code{"AIC"}, \code{"cAIC"},
#' \code{"mAIC"}, \code{"BIC"}, \code{"mBIC"}, \code{"EBIC"}, and
#' \code{"logLik"}. This option has no effect for objects fit using
#' \code{cv_curegmifs} or \code{cv_cureem}.
#' @param xlab  a default x-axis label will be used which can be changed by
#' specifying a user-defined x-axis label.
#' @param ylab a default y-axis label will be used which can be changed by
#' specifying a user-defined y-axis label.
#' @param main a default main title will be used which can be changed by
#' specifying a user-defined main title. This option is not used for
#' \code{cv_curegmifs} or \code{cv_cureem} fitted objects.
#' @param ... other arguments.
#'
#' @return this function has no returned value but is called for its side
#' effects
#'
#' @seealso \code{\link{curegmifs}}, \code{\link{cureem}},
#' \code{\link{coef.mixturecure}}, \code{\link{summary.mixturecure}},
#' \code{\link{predict.mixturecure}}
#' @keywords methods
#' @method plot mixturecure
#' @import ggplot2
#' @import graphics
#' @importFrom ggpubr ggarrange
#' @export
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'                    data = training, x_latency = training,
#'                    model = "weibull", thresh = 1e-4, maxit = 2000,
#'                    epsilon = 0.01, verbose = FALSE)
#' plot(fit)
plot.mixturecure <-
  function(x, type = c("trace", "AIC", "BIC", "logLik", "cAIC", "mAIC",
                       "mBIC", "EBIC"), xlab = NULL, ylab = NULL, main = NULL,
           ...) {
    #type <- c("trace", "AIC", "BIC", "logLik", "cAIC", "mAIC", "mBIC", "EBIC")[pmatch(type,
    #                                                                                  c("trace", "AIC", "BIC", "logLik", "cAIC", "mAIC", "mBIC", "EBIC"))]
    #if (!(type%in%c("trace", "AIC", "BIC", "logLik", "cAIC", "mAIC", "mBIC", "EBIC")))
    #  stop("type must be either 'trace', 'AIC', 'BIC', 'logLik', 'cAIC', 'mAIC', 'mBIC', or 'EBIC' ")
    type <- match.arg(type)
    if (!x$cv) {
      p <- dim(x$x.incidence)[2] + dim(x$x_latency)[2]
      if (type %in% c("AIC", "BIC", "cAIC", "mAIC", "mBIC", "EBIC")) {
        if (!is.null(x$x.incidence)) {
          vars_inc <- apply(x$b_path, 1, function(x) sum(x != 0))
        } else {
          vars_inc <- 0
        }
        if (!is.null(x$x_latency)) {
          vars_lat <- apply(x$beta_path, 1, function(x) sum(x != 0))
        } else {
          vars_lat <- 0
        }
        if (x$model == "weibull") {
          df <- vars_inc + vars_lat + 3
        } else if (x$model == "exponential") {
          df <- vars_inc + vars_lat + 2
        } else if (x$model == "cox") {
          df <- vars_inc + vars_lat + 1
        }
      }
      if (is.null(xlab))
        xlab <- "Step"
      if (is.null(ylab)) {
        if (type == "AIC") {
          ylab <- "AIC"
        } else if (type == "cAIC") {
          ylab <- "cAIC"
        } else if (type == "mAIC") {
          ylab <- "mAIC"
        } else if (type == "BIC") {
          ylab <- "BIC"
        } else if (type == "mBIC") {
          ylab <- "mBIC"
        } else if (type == "EBIC") {
          ylab <- "EBIC"
        } else if (type == "trace") {
          ylab <- expression(hat(beta))
        } else if (type == "logLik") {
          ylab <- "logLikelihood"
        }
      }
      if (x$method == "EM")
        logLik <- x$logLik_inc + x$logLik_lat
      else
        logLik <- x$logLik
      if (!is.null(x$x_latency) && !is.null(x$x.incidence)) {
        coef <- cbind(b = x$b_path, beta = x$beta_path)
      } else if (is.null(x$x_latency) && !is.null(x$x.incidence)) {
        coef <- x$b_path
      } else if (!is.null(x$x_latency) && is.null(x$x.incidence)) {
        coef <- x$beta_path
      }
      if (type == "trace") {
        graphics::matplot(coef, ylab = ylab, xlab = xlab, type = "l")
      } else if (type %in% c("AIC", "cAIC", "mAIC")) {
        AIC <- 2 * df - 2 * logLik
        cAIC <- AIC + (2 * df * (df + 1)) / (length(x$y) - df - 1)
        mAIC <- (2 + 2 * log(p / .5)) * df - 2 * logLik
        if (type == "AIC") {
          plot(AIC, xlab = xlab, ylab = ylab)
        } else if (type == "cAIC") {
          plot(cAIC, xlab = xlab, ylab = ylab)
        } else {
          plot(mAIC, xlab = xlab, ylab = ylab)
        }
      } else if (type %in% c("BIC", "mBIC", "EBIC")) {
        BIC <- df * (log(length(x$y))) - 2 * logLik
        mBIC <- df * (log(length(x$y)) + 2 * log(p / 4)) - 2 * logLik
        EBIC <- log(length(x$y)) * df + 2 * (1 - .5) * log(choose(p, df)) -
          2 * logLik
        if (type == "BIC") {
          plot(BIC, xlab = xlab, ylab = ylab)
        } else if (type == "EBIC") {
          plot(EBIC, xlab = xlab, ylab = ylab)
        } else if (type == "mBIC") {
          plot(mBIC, xlab = xlab, ylab = ylab)
        }
      } else if (type == "logLik") {
        plot(logLik, xlab = xlab, ylab = ylab)
      }
      if (is.null(main)) {
        graphics::title(paste(type, "from MCM ", x$method, " fit",
                              sep = " "))
      } else {
        graphics::title(main)
      }
    } else {
      if (is.null(xlab)) xlab <- "Predictors"
      if (is.null(ylab)) ylab <- "Coefficient Estimate"
      names(x$b) <- colnames(x$x.incidence)
      names(x$beta) <- colnames(x$x_latency)
      if (!is.null(x$x.incidence)) {
        if (sum(x$b > 0) > 0) {
          b <- x$b[x$b != 0]
        } else {
          b <- x$b
        }
      }
      if (!is.null(x$x_latency)) {
        if (sum(x$beta > 0) > 0) {
          beta <- x$beta[x$beta != 0]
        } else {
          beta <- x$beta
        }
      }
      y <- NULL
      if (!is.null(x$x.incidence))
        data_b <- data.frame(x = seq_along(b), y = b, vars = names(b),
                             model_component = "Incidence")
      if (!is.null(x$x_latency))
        data_beta <- data.frame(x = seq_along(beta), y = beta,
                                vars = names(beta), model_component = "Latency")
      if (!is.null(x$x.incidence) && !is.null(x$x_latency)) {
        bplot <- ggplot(data_b, aes(x = vars, y = y)) + geom_point() +
          geom_segment(aes(x = vars, xend = vars, y = 0, yend = y)) +
          xlab(xlab) + ylab(ylab) + ggtitle("Incidence") +
          theme(axis.text.x = element_text(angle = 45))
        betaplot <- ggplot(data_beta, aes(x = vars, y = y)) + geom_point() +
          geom_segment(aes(x = vars, xend = vars, y = 0, yend = y)) +
          xlab(xlab) + ylab(ylab) + ggtitle("Latency") +
          theme(axis.text.x = element_text(angle = 45))
        ggpubr::ggarrange(bplot, betaplot, nrow = 2)
      } else if (is.null(x$x.incidence) && !is.null(x$x_latency)) {
        ggplot(data_beta, aes(x = vars, y = y)) + geom_point() +
          geom_segment(aes(x = vars, xend = vars, y = 0, yend = y)) +
          xlab(xlab) + ylab(ylab) + ggtitle("Latency") +
          theme(axis.text.x = element_text(angle = 45))
      } else if (!is.null(x$x.incidence) && is.null(x$x_latency)) {
        ggplot(data_b, aes(x = vars, y = y)) + geom_point() +
          geom_segment(aes(x = vars, xend = vars, y = 0, yend = y)) +
          xlab(xlab) + ylab(ylab) + ggtitle("Incidence") +
          theme(axis.text.x = element_text(angle = 45))
      }
    }
  }
