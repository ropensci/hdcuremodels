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
#' @param type a case-sensitive parameter indicating what to plot on the y-axis.
#' The complete list of options are:
#' \itemize{
#'     \item \code{"trace"} plots the coefficient path for the fitted object
#'     (default).
#'     \item \code{"AIC"} plots the AIC against step of model fit.
#'     \item \code{"mAIC"} plots the modified AIC against step of model fit.
#'     \item \code{"cAIC"} plots the corrected AIC against step of model fit.
#'     \item \code{"BIC"}, plots the BIC against step of model fit.
#'     \item \code{"mBIC"} plots the modified BIC against step of model fit.
#'     \item \code{"EBIC"} plots the extended BIC against step of model fit.
#'     \item \code{"logLik"} plots the log-likelihood against step of model fit.
#'   }
#' This option has no effect for objects fit using
#' \code{cv_curegmifs} or \code{cv_cureem}.
#' @param xlab  a default x-axis label will be used which can be changed by
#' specifying a user-defined x-axis label.
#' @param ylab a default y-axis label will be used which can be changed by
#' specifying a user-defined y-axis label.
#' @param label logical. If TRUE the variable names will appear in a legend.
#' Applicable only when \code{type = "trace"}. Be reminded that this works well
#' only for small to moderate numbers of variables. For many predictors, the
#' plot will be cluttered. The variables may be more easily identified using
#' the \code{coef} function indicating the step of interest.
#' @param main a default main title will be used which can be changed by
#' specifying a user-defined main title. This option is not used for
#' \code{cv_curegmifs} or \code{cv_cureem} fitted objects.
#' @param ... other arguments.
#'
#' @return this function has no returned value but is called for its side
#' effects
#'
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {RE6.0} *Model objects returned by Regression Software (see* **RE4***) should have default `plot` methods, either through explicit implementation, extension of methods for existing model objects, or through ensuring default methods work appropriately.*
#' @srrstats {RE6.1} *Where the default `plot` method is **NOT** a generic `plot` method dispatched on the class of return objects (that is, through an S3-type `plot.<myclass>` function or equivalent), that method dispatch (or equivalent) should nevertheless exist in order to explicitly direct users to the appropriate function.*
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
#' withr::local_seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'   data = training, x_latency = training,
#'   model = "weibull", thresh = 1e-4, maxit = 2000,
#'   epsilon = 0.01, verbose = FALSE
#' )
#' plot(fit)
plot.mixturecure <-
  function(x, type = c(
             "trace", "AIC", "BIC", "logLik", "cAIC", "mAIC",
             "mBIC", "EBIC"), xlab = NULL, ylab = NULL,
              label = FALSE, main = NULL, ...) {
    type <- match.arg(type)
    if (!x$cv) {
      if (is.null(xlab)) {
        xlab <- "Step"
      }
      if (type == "trace") {
        if (is.null(ylab)) {
          ylab <- expression(hat(beta))
        }
        if (!is.null(x$x_latency) && !is.null(x$x_incidence)) {
          colnames(x$b_path) <- paste0("I_", colnames(x$b_path))
          colnames(x$beta_path) <- paste0("L_", colnames(x$beta_path))
          coef <- cbind(b = x$b_path, beta = x$beta_path)
        } else if (is.null(x$x_latency) && !is.null(x$x_incidence)) {
          colnames(x$b_path) <- paste0("I_", colnames(x$b_path))
          coef <- x$b_path
        } else if (!is.null(x$x_latency) && is.null(x$x_incidence)) {
          colnames(x$beta_path) <- paste0("L_", colnames(x$beta_path))
          coef <- x$beta_path
        }
        graphics::matplot(coef, ylab = ylab, xlab = xlab, type = "l",
                    lty = rep(c(1,2), c(dim(x$b_pat)[2], dim(x$beta_path)[2])),
                    col = rep(1:(dim(x$b_path)[2]+dim(x$beta_path)[2])))
        if (label == TRUE)
          legend("topright", legend = colnames(coef),
                 lty = rep(c(1,2), c(dim(x$b_pat)[2], dim(x$beta_path)[2])),
                 cex = 0.6, col = rep(1:(dim(x$b_path)[2]+dim(x$beta_path)[2])))
      } else {
        select <- select_model(x, type)
        if (type == "AIC") {
        if (is.null(ylab)) {
          ylab <- "AIC"
        }
        plot(select$AIC, xlab = xlab, ylab = ylab)
      } else if (type == "cAIC") {
        if (is.null(ylab)) {
          ylab <- "cAIC"
        }
        plot(select$cAIC, xlab = xlab, ylab = ylab)
      } else if (type == "mAIC") {
        if (is.null(ylab)) {
          ylab <- "mAIC"
        }
        plot(select$mAIC, xlab = xlab, ylab = ylab)
      } else if (type == "BIC") {
        if (is.null(ylab)) {
          ylab <- "BIC"
        }
        plot(select$BIC, xlab = xlab, ylab = ylab)
      } else if (type == "mBIC") {
        if (is.null(ylab)) {
          ylab <- "mBIC"
        }
        plot(select$mBIC, xlab = xlab, ylab = ylab)
      } else if (type == "EBIC") {
        if (is.null(ylab)) {
          ylab <- "EBIC"
        }
        plot(select$EBIC, xlab = xlab, ylab = ylab)
      } else if (type == "logLik") {
        if (is.null(ylab)) {
          ylab <- "logLikelihood"
        }
        plot(select$logLik, xlab = xlab, ylab = ylab)
      }
      if (is.null(main)) {
        graphics::title(paste(type, "from MCM ", x$method, " fit",
          sep = " "
        ))
      } else {
        graphics::title(main)
      }
    }
    } else {
      if (is.null(xlab)) xlab <- "Predictors"
      if (is.null(ylab)) ylab <- "Coefficient Estimates"
      names(x$b) <- colnames(x$x_incidence)
      names(x$beta) <- colnames(x$x_latency)
      if (!is.null(x$x_incidence)) {
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
      if (!is.null(x$x_incidence)) {
        data_b <- data.frame(
          x = seq_along(b), y = b, vars = names(b),
          model_component = "Incidence"
        )
      }
      if (!is.null(x$x_latency)) {
        data_beta <- data.frame(
          x = seq_along(beta), y = beta,
          vars = names(beta), model_component = "Latency"
        )
      }
      if (!is.null(x$x_incidence) && !is.null(x$x_latency)) {
        bplot <- ggplot(data_b, aes(x = vars, y = y)) +
          geom_point() +
          geom_segment(aes(x = vars, xend = vars, y = 0, yend = y)) +
          xlab(xlab) +
          ylab(ylab) +
          ggtitle("Incidence") +
          theme(axis.text.x = element_text(angle = 45))
        betaplot <- ggplot(data_beta, aes(x = vars, y = y)) +
          geom_point() +
          geom_segment(aes(x = vars, xend = vars, y = 0, yend = y)) +
          xlab(xlab) +
          ylab(ylab) +
          ggtitle("Latency") +
          theme(axis.text.x = element_text(angle = 45))
        ggpubr::ggarrange(bplot, betaplot, nrow = 2)
      } else if (is.null(x$x_incidence) && !is.null(x$x_latency)) {
        ggplot(data_beta, aes(x = vars, y = y)) +
          geom_point() +
          geom_segment(aes(x = vars, xend = vars, y = 0, yend = y)) +
          xlab(xlab) +
          ylab(ylab) +
          ggtitle("Latency") +
          theme(axis.text.x = element_text(angle = 45))
      } else if (!is.null(x$x_incidence) && is.null(x$x_latency)) {
        ggplot(data_b, aes(x = vars, y = y)) +
          geom_point() +
          geom_segment(aes(x = vars, xend = vars, y = 0, yend = y)) +
          xlab(xlab) +
          ylab(ylab) +
          ggtitle("Incidence") +
          theme(axis.text.x = element_text(angle = 45))
      }
    }
  }
