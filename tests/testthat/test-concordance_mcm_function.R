test_that("concordance_mcm works correctly", {
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
  expect_equal(round(concordance_mcm(fit, newdata = testing,
                                     model_select = "cAIC"), 7), 0.7455192)
  expect_equal(round(concordance_mcm(fit, cure_cutoff = 3,
                                     model_select = "cAIC"), 6), 0.8598)
  expect_equal(round(concordance_mcm(fit, model_select = "cAIC"), 7), 0.8627683)
  expect_equal(round(concordance_mcm(fit, model_select = "AIC"), 7), 0.8627683)
  expect_equal(round(concordance_mcm(fit, model_select = "mAIC"), 6), 0.546498)
  expect_equal(round(concordance_mcm(fit, model_select = "BIC"), 7), 0.8627683)
  expect_equal(round(concordance_mcm(fit, model_select = "mBIC"), 6), 0.546498)
  expect_equal(round(concordance_mcm(fit, model_select = "EBIC"), 7), 0.8614357)
  expect_equal(round(concordance_mcm(fit, model_select = "logLik"), 7), 0.8627683)
  expect_error(concordance_mcm(fit, model_select = "aic"))
  expect_error(concordance_mcm(fit, model_select = "bic"))
  expect_error(concordance_mcm("x"))
  expect_error(concordance_mcm(1))

})
