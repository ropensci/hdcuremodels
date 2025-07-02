#' @srrstats {G2.0} *Implement assertions on lengths of inputs, particularly through asserting that inputs expected to be single- or multi-valued are indeed so.*
#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*

test_that("multiplication works", {
  library(survival)
  withr::local_seed(1234)
  temp <- generate_cure_data(n = 200, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  training$group <- gl(2, 75)
  km_fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
  output <- cure_estimate(km_fit)
  output %>% expect_length(1)
  expect_equal(round(output, 7), 0.4866733)
  output %>% expect_type("double")
  expect_error(cure_estimate(survdiff(Surv(Time, Censor) ~ 1, data = training)))
  expect_visible(cure_estimate(km_fit))

  km_fit <- survfit(Surv(Time, Censor) ~ group, data = training)
  output <- cure_estimate(km_fit)
  output %>% expect_length(3)
  expect_equal(round(output[1, 2], 7), 0.5357901)
  output %>% expect_type("list")
  expect_error(cure_estimate(survdiff(Surv(Time, Censor) ~ group, data = training)))
  expect_visible(cure_estimate(km_fit))

  fit.lm <- lm(Time ~ Censor, data = training)
  expect_error(cure_estimate(fit.lm), "Error: object must be a survfit object")

})
