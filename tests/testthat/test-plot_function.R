test_that("plot.mixturecure function works correctly", {
  withr::local_seed(1234)
  temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  fit <- curegmifs(Surv(Time, Censor) ~ .,
                   data = training, x_latency = training,
                   model = "weibull", thresh = 1e-4, maxit = 2000,
                   epsilon = 0.01, verbose = FALSE
  )
  expect_invisible(plot(fit, xlab = "Step", ylab = "Log-likelihood",
                main = "trace"))
  expect_invisible(plot(fit, type = "AIC"))
  expect_invisible(plot(fit, type = "BIC"))
  expect_invisible(plot(fit, type = "logLik"))
  expect_invisible(plot(fit, type = "mAIC"))
  expect_invisible(plot(fit, type = "mBIC"))
  expect_invisible(plot(fit, type = "EBIC"))

  withr::local_seed(123)
  temp <- generate_cure_data(n = 100, j = 15, n_true = 3, a = 1.8, rho = 0.2)
  training <- temp$training
  fit.cv <- cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training,
                         x_latency = training, fdr_control = FALSE,
                         maxit = 450, epsilon = 0.01, n_folds = 2,
                         seed = 23, verbose = TRUE
  )
  expect_visible(plot(fit.cv))

  fit.cv.fdr <- cv_cureem(Surv(Time, Censor) ~ .,
                          data = training,
                          x_latency = training, model = "weibull", penalty = "lasso",
                          fdr_control = TRUE, grid_tuning = FALSE, nlambda_inc = 10,
                          nlambda_lat = 10, n_folds = 2, seed = 23, verbose = TRUE
  )
  expect_silent(plot(fit.cv.fdr))
  fit.em <- cureem(Surv(Time, Censor) ~ .,
                data = training, x_latency = training,
                model = "cox", penalty = "lasso", lambda_inc = 0.1,
                lambda_lat = 0.1, gamma_inc = 6, gamma_lat = 10
  )
  expect_silent(plot(fit.em))
})
