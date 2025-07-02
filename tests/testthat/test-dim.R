test_that("dim function works correctly", {
  library(survival)
  withr::local_seed(1234)
  temp <- generate_cure_data(n = 80, j = 100, n_true = 10, a = 1.8)
  training <- temp$training
  fit <- cureem(Surv(Time, Censor) ~ .,
                data = training, x_latency = training,
                model = "cox", penalty = "lasso", lambda_inc = 0.1,
                lambda_lat = 0.1, gamma_inc = 6, gamma_lat = 10
  )
  expect_equal(dim(fit)[["nobs"]], 60)
  expect_equal(dim(fit)[["p_incidence"]], 102)
  expect_equal(dim(fit)[["p_latency"]], 102)
})
