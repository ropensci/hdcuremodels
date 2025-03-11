test_that("print works correctly", {
  library(survival)
  set.seed(1234)
  temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  fit <- curegmifs(Surv(Time, Censor) ~ .,
                   data = training, x_latency = training,
                   model = "weibull", thresh = 1e-4, maxit = 2000,
                   epsilon = 0.01, verbose = FALSE
  )
  result <- print(fit)
  expect_true(is.list(result))
  result %>% expect_length(15)
  expect_equal(names(result), c("b_path", "beta_path", "b0_path", "rate_path",
                             "logLik", "x_incidence", "x_latency", "y",
                             "model", "scale", "method", "call",
                             "alpha_path", "cv", "warning" ))
  expect_equal(round(mean(result$alpha_path), 6), 1.132515)

  expect_warning(curegmifs(Surv(Time, Censor) ~ .,
                   data = training, x_latency = training,
                   model = "exponential", thresh = 1e-4, maxit = 2000,
                   epsilon = 0.01, verbose = FALSE
  ))

  fit.cv <- cv_cureem(Surv(Time, Censor) ~ .,
                      data = training,
                      x_latency = training, fdr_control = FALSE,
                      grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
                      n_folds = 2, seed = 23, verbose = TRUE
  )
  result <- print(fit.cv)
  expect_true(is.list(result))
  result %>% expect_length(18)
  expect_visible(result)

  fit.cv.fdr <- cv_cureem(Surv(Time, Censor) ~ .,
                          data = training,
                          x_latency = training, model = "weibull", penalty = "lasso",
                          fdr_control = TRUE, grid_tuning = FALSE, nlambda_inc = 10,
                          nlambda_lat = 10, n_folds = 2, seed = 23, verbose = TRUE
  )
  result <- print(fit.cv.fdr)
  expect_true(is.list(result))
  result %>% expect_length(17)
  expect_visible(result)

})
