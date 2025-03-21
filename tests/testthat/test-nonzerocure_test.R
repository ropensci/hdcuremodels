#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstats {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*

test_that("nonzerocure_test works", {
  library(survival)
  withr::local_seed(1234)
  temp <- generate_cure_data(n = 200, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  training$group <- gl(2, 75)
  km_fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
  output <- nonzerocure_test(km_fit)
  output %>% expect_type("list")
  output %>% expect_length(4)
  expect_equal(round(output$p_value, 3), 0.039)

  fit.lm <- lm(Time ~ Censor, data = training)
  expect_error(nonzerocure_test(fit.lm), "Error: object must be a survfit object")
})
