#' Print the contents of a mixture cure fitted object
#'
#' @description
#' This function prints the names of the list objects from a \code{curegmifs},
#' \code{cureem}, \code{cv_cureem}, or \code{cv_curegmifs} fitted model.
#'
#' @param x a \code{mixturecure} object resulting from \code{curegmifs},
#' \code{cureem}, \code{cv_cureem}, or \code{cv_curegmifs}.
#' @param ... other arguments.
#'
#' @note The contents of an \code{mixturecure} fitted object differ depending
#' upon whether the EM (\code{cureem}) or GMIFS (\code{curegmifs}) algorithm is
#' used for model fitting. Also, the output differs depending upon whether
#' \code{x_latency} is specified in the model (i.e., variables are included in
#' the latency portion of the model fit) or only \code{terms} on the right hand
#' side of the equation are included (i.e., variables are included in the
#' incidence portion of the model).
#'
#' @seealso \code{\link{curegmifs}}, \code{\link{cureem}},
#' \code{\link{coef.mixturecure}}, \code{\link{summary.mixturecure}},
#' \code{\link{plot.mixturecure}}, \code{\link{predict.mixturecure}}
#'
#' @return names of the objects in a mixturecure object fit using \code{cureem},
#' \code{curegmifs}, \code{cv_cureem}, or \code{cv_curegmifs}.
#' @export
#'
#' @keywords methods
#' @method print mixturecure
#' @examples
#' library(survival)
#' set.seed(1234)
#' temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
#' training <- temp$training
#' fit <- curegmifs(Surv(Time, Censor) ~ .,
#'                        data = training, x_latency = training,
#'                        model = "weibull", thresh = 1e-4, maxit = 2000,
#'                        epsilon = 0.01, verbose = FALSE)
#' print(fit)
print.mixturecure <- function(x, ...) {
  print(names(x))
}
