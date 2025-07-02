#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*

test_that("auc_mcm works correctly", {
  library(survival)
  withr::local_seed(1234)
  temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  testing <- temp$testing
  fit <- curegmifs(Surv(Time, Censor) ~ .,
                   data = training, x_latency = training,
                   model = "weibull", thresh = 1e-4, maxit = 2000,
                   epsilon = 0.01, verbose = FALSE
  )
  expect_equal(round(auc_mcm(fit, cure_cutoff = 3, model_select = "cAIC"), 7),
               0.7762012)
  expect_equal(round(auc_mcm(fit, newdata = testing, cure_cutoff = 3,
                             model_select = "cAIC"), 7), 0.6709119)
  expect_equal(round(auc_mcm(fit, model_select = "cAIC"), 7), 0.8141595)
  expect_equal(round(auc_mcm(fit, model_select = "AIC"), 7), 0.8263137)
  expect_equal(round(auc_mcm(fit, model_select = "mAIC"), 6), 0.678221)
  expect_equal(round(auc_mcm(fit, model_select = "BIC"), 7), 0.8141595)
  expect_equal(round(auc_mcm(fit, model_select = "mBIC"), 6), 0.678221)
  expect_equal(round(auc_mcm(fit, model_select = "EBIC"), 6), 0.813883)
  expect_equal(round(auc_mcm(fit, model_select = "logLik"), 7), 0.8285604)
  expect_equal(round(auc_mcm(fit, model_select = 22), 7), 0.691254)
  expect_error(auc_mcm(fit, model_select = "aic"))
  expect_error(auc_mcm(fit, model_select = "bic"))
  expect_error(auc_mcm("x"))
  expect_error(auc_mcm(1))

  fit.lm <- lm(Time ~ Censor, data = training)
  expect_error(auc_mcm(fit.lm), "Error: class of object must be mixturecure")
})
