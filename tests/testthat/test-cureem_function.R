test_that("cureem works correctly", {
# validate function output
  library(survival)
  set.seed(1234)
  temp <- generate_cure_data(n = 80, j = 100, n_true = 10, a = 1.8)
  training <- temp$training
  fit <- cureem(Surv(Time, Censor) ~ .,
                data = training, x_latency = training,
                model = "cox", penalty = "lasso", lambda_inc = 0.1,
                lambda_lat = 0.1, gamma_inc = 6, gamma_lat = 10
  )
  expect_equal(names(fit), c("b_path", "beta_path", "b0_path", "logLik_inc",
                             "logLik_lat", "x_incidence", "x_latency", "y",
                             "model", "scale", "method", "call", "cv" ))
  expect_true(class(fit) == "mixturecure")
  expect_true(class(fit$b_path)[1] == "matrix")
  expect_true(class(fit$beta_path)[1] == "matrix")
  expect_true(class(fit$b0_path) == "numeric")
  expect_true(class(fit$logLik_inc) == "numeric")
  expect_true(class(fit$logLik_lat) == "numeric")


  fit <- cureem(Surv(Time, Censor) ~ .,
                data = training, x_latency = training,
                model = "cox", penalty = "MCP", lambda_inc = 0.1,
                lambda_lat = 0.1, gamma_inc = 6, gamma_lat = 10
  )
  expect_equal(names(fit), c("b_path", "beta_path", "b0_path", "logLik_inc",
                             "logLik_lat", "x_incidence", "x_latency", "y",
                             "model", "scale", "method", "call", "cv" ))
  expect_true(class(fit) == "mixturecure")
  expect_true(class(fit$b_path)[1] == "matrix")
  expect_true(class(fit$beta_path)[1] == "matrix")
  expect_true(class(fit$b0_path) == "numeric")
  expect_true(class(fit$logLik_inc) == "numeric")
  expect_true(class(fit$logLik_lat) == "numeric")

  fit <- cureem(Surv(Time, Censor) ~ .,
                data = training, x_latency = training,
                model = "cox", penalty = "SCAD", lambda_inc = 0.1,
                lambda_lat = 0.1, gamma_inc = 6, gamma_lat = 10
  )
  expect_equal(names(fit), c("b_path", "beta_path", "b0_path", "logLik_inc",
                             "logLik_lat", "x_incidence", "x_latency", "y",
                             "model", "scale", "method", "call", "cv" ))
  expect_true(class(fit) == "mixturecure")
  expect_true(class(fit$b_path)[1] == "matrix")
  expect_true(class(fit$beta_path)[1] == "matrix")
  expect_true(class(fit$b0_path) == "numeric")
  expect_true(class(fit$logLik_inc) == "numeric")
  expect_true(class(fit$logLik_lat) == "numeric")


  fit <- cureem(Surv(Time, Censor) ~ .,
                data = training, x_latency = training,
                model = "weibull", penalty = "lasso", lambda_inc = 0.1,
                lambda_lat = 0.1, gamma_inc = 6, gamma_lat = 10
  )
  expect_equal(names(fit), c("b_path", "beta_path", "b0_path", "logLik_inc",
                             "logLik_lat", "x_incidence", "x_latency", "y",
                             "model", "scale", "method", "call", "rate",
                             "alpha", "cv" ))
  expect_true(class(fit) == "mixturecure")
  expect_true(class(fit$b_path)[1] == "matrix")
  expect_true(class(fit$beta_path)[1] == "matrix")
  expect_true(class(fit$b0_path) == "numeric")
  expect_true(class(fit$logLik_inc) == "numeric")
  expect_true(class(fit$logLik_lat) == "numeric")

  fit <- cureem(Surv(Time, Censor) ~ .,
                data = training, x_latency = training,
                model = "exponential", penalty = "lasso", lambda_inc = 0.1,
                lambda_lat = 0.1, gamma_inc = 6, gamma_lat = 10
  )
  expect_equal(names(fit), c("b_path", "beta_path", "b0_path", "logLik_inc",
                             "logLik_lat", "x_incidence", "x_latency", "y",
                             "model", "scale", "method", "call", "rate",
                              "cv" ))
  expect_true(class(fit) == "mixturecure")
  expect_true(class(fit$b_path)[1] == "matrix")
  expect_true(class(fit$beta_path)[1] == "matrix")
  expect_true(class(fit$b0_path) == "numeric")
  expect_true(class(fit$logLik_inc) == "numeric")
  expect_true(class(fit$logLik_lat) == "numeric")

})
