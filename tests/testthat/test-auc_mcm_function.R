test_that("auc_mcm works correctly", {
  library(survival)
  set.seed(1234)
  temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  testing <- temp$testing
  fit <- curegmifs(Surv(Time, Censor) ~ .,
                   data = training, x_latency = training,
                   model = "weibull", thresh = 1e-4, maxit = 2000,
                   epsilon = 0.01, verbose = FALSE
  )
  expect_equal(round(auc_mcm(fit, cure_cutoff = 3, model_select = "cAIC"), 7),
               0.8403799)
  expect_equal(round(auc_mcm(fit, newdata = testing, cure_cutoff = 3,
                             model_select = "cAIC"), 7), 0.8066537)
  expect_equal(round(auc_mcm(fit, model_select = "cAIC"), 7), 0.8754937)
  expect_equal(round(auc_mcm(fit, model_select = "AIC"), 7), 0.8754937)
  expect_equal(round(auc_mcm(fit, model_select = "mAIC"), 6), 0.706476)
  expect_equal(round(auc_mcm(fit, model_select = "BIC"), 7), 0.8754937)
  expect_equal(round(auc_mcm(fit, model_select = "mBIC"), 6), 0.706476)
  expect_equal(round(auc_mcm(fit, model_select = "EBIC"), 6), 0.872958)
  expect_equal(round(auc_mcm(fit, model_select = "logLik"), 7), 0.8754937)
  expect_equal(round(auc_mcm(fit, model_select = 22), 7), 0.7783467)
  expect_error(auc_mcm(fit, model_select = "aic"))
  expect_error(auc_mcm(fit, model_select = "bic"))
  expect_error(auc_mcm("x"))
  expect_error(auc_mcm(1))

  fit.lm <- lm(Time ~ Censor, data = training)
  expect_error(auc_mcm(fit.lm), "Error: class of object must be mixturecure")
})
