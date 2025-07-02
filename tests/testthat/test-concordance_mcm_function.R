#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*

test_that("concordance_mcm works correctly", {
  library(survival)
  withr::local_seed(234)
  temp <- generate_cure_data(n = 100, j = 10, n_true = 8, a = 2)
  training <- temp$training
  testing <- temp$testing
  fit <- curegmifs(Surv(Time, Censor) ~ .,
                   data = training, x_latency = training,
                   model = "weibull", thresh = 1e-4, maxit = 2000,
                   epsilon = 0.01, verbose = FALSE
  )
  expect_equal(round(concordance_mcm(fit, newdata = testing,
                                     model_select = "cAIC"), 7), 0.7320809)
  expect_equal(round(concordance_mcm(fit, cure_cutoff = 3,
                                     model_select = "cAIC"), 6), 0.883804)
  expect_equal(round(concordance_mcm(fit, model_select = "cAIC"), 7), 0.8830314)
  expect_equal(round(concordance_mcm(fit, model_select = "AIC"), 7), 0.8830314)
  expect_equal(round(concordance_mcm(fit, model_select = "mAIC"), 6), 0.662166)
  expect_equal(round(concordance_mcm(fit, model_select = "BIC"), 7), 0.8830314)
  expect_equal(round(concordance_mcm(fit, model_select = "mBIC"), 6), 0.662166)
  expect_equal(round(concordance_mcm(fit, model_select = "EBIC"), 7), 0.8581016)
  expect_equal(round(concordance_mcm(fit, model_select = "logLik"), 7), 0.8830314)
  expect_error(concordance_mcm(fit, model_select = "aic"))
  expect_error(concordance_mcm(fit, model_select = "bic"))
  expect_error(concordance_mcm("x"))
  expect_error(concordance_mcm(1))

  fit.lm <- lm(Time ~ Censor, data = training)
  expect_error(concordance_mcm(fit.lm), "Error: class of object must be mixturecure")

})
