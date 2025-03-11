test_that("summary function works correctly", {
  library(survival)
  set.seed(1234)
  temp <- generate_cure_data(n = 80, j = 100, n_true = 10, a = 1.8)
  training <- temp$training
  fit <- cureem(Surv(Time, Censor) ~ .,
                data = training, x_latency = training,
                model = "cox", penalty = "lasso", lambda_inc = 0.1,
                lambda_lat = 0.1, gamma_inc = 6, gamma_lat = 10
  )
  expect_invisible(summary(fit))
  expect_null(summary(fit))
  fit.cv <- cv_cureem(Surv(Time, Censor) ~ .,
                      data = training,
                      x_latency = training, fdr_control = FALSE,
                      grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
                      n_folds = 2, seed = 23, verbose = TRUE
  )
  expect_invisible(summary(fit.cv))
  expect_null(summary(fit.cv))

  fit.lm <- lm(Time ~ Censor, data = training)
  expect_error(summary.mixturecure(fit.lm))
})
