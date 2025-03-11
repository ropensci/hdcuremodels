#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*

test_that("cv_curegmifs function works correctly", {
  library(survival)
  set.seed(123)
  temp <- generate_cure_data(n = 100, j = 15, n_true = 3, a = 1.8, rho = 0.2)
  training <- temp$training

  fit.cv <- cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training,
                         x_latency = training, fdr_control = FALSE,
                         maxit = 450, epsilon = 0.01, n_folds = 2,
                         seed = 23, verbose = TRUE
  )
  expect_error(cv_curegmifs(Time ~ ., data = training, x_latency = training))
  expect_error(cv_curegmifs(Censor ~ ., data = training, x_latency = training))

  fit.cv %>% expect_s3_class("mixturecure")
  expect_setequal(names(fit.cv), c("selected_step_inc", "selected_step_lat",
                                   "b0", "b", "beta", "rate", "logLik", "max_c",
                                   "x_incidence", "x_latency", "y", "model",
                                   "scale", "method", "alpha", "cv", "call",
                                   "fdr_control"))
  fit.cv$selected_step_inc %>% expect_type("integer")
  fit.cv$selected_step_lat %>% expect_type("integer")
  fit.cv$b0 %>% expect_type("double")
  fit.cv$b %>% expect_type("double")
  fit.cv$beta %>% expect_type("double")
  fit.cv$rate %>% expect_type("double")
  fit.cv$logLik %>% expect_type("double")
  fit.cv$max_c %>% expect_type("double")
  fit.cv$x_incidence %>% expect_type("double")
  fit.cv$x_latency %>% expect_type("double")
  fit.cv$y %>% expect_type("double")
  fit.cv$model %>% expect_type("character")
  fit.cv$scale %>% expect_type("logical")
  fit.cv$method %>% expect_type("character")
  fit.cv$alpha %>% expect_type("double")
  fit.cv$cv %>% expect_type("logical")
  fit.cv$fdr_control %>% expect_type("logical")
  fit.cv$call %>% expect_type("language")

  fit.cv$b %>% expect_length(dim(fit.cv$x_incidence)[2])
  fit.cv$beta %>% expect_length(dim(fit.cv$x_latency)[2])
  fit.cv$b0 %>% expect_length(1)
  fit.cv$alpha %>% expect_length(1)
  fit.cv$rate %>% expect_length(1)
  fit.cv$selected_step_inc %>% expect_length(1)
  fit.cv$selected_step_lat %>% expect_length(1)
  fit.cv$logLik %>% expect_length(1)
  fit.cv$max_c %>% expect_length(1)
  expect_equal(round(fit.cv$rate, 6), 2.502845)
  expect_equal(round(fit.cv$alpha, 7), 0.7003362)

  fit.cv <- cv_curegmifs(Surv(Time, Censor) ~ ., data = training,
                         penalty_factor_inc = rep(c(0, 1), c(1, 11)),
                         measure_inc = "auc",
                         x_latency = training, fdr_control = FALSE,
                         maxit = 450, epsilon = 0.01, n_folds = 2,
                         seed = 23, verbose = FALSE, parallel = FALSE
  )
  fit.cv %>% expect_s3_class("mixturecure")
  expect_setequal(names(fit.cv), c("selected_step_inc", "selected_step_lat",
                                   "b0", "b", "beta", "rate", "logLik", "max_c",
                                   "x_incidence", "x_latency", "y", "model",
                                   "scale", "method", "alpha", "max.auc", "cv",
                                   "call", "fdr_control"))
  fit.cv$selected_step_inc %>% expect_type("integer")
  fit.cv$selected_step_lat %>% expect_type("integer")
  fit.cv$b0 %>% expect_type("double")
  fit.cv$b %>% expect_type("double")
  fit.cv$beta %>% expect_type("double")
  fit.cv$rate %>% expect_type("double")
  fit.cv$logLik %>% expect_type("double")
  fit.cv$max_c %>% expect_type("double")
  fit.cv$max.auc %>% expect_type("double")
  fit.cv$x_incidence %>% expect_type("double")
  fit.cv$x_latency %>% expect_type("double")
  fit.cv$y %>% expect_type("double")
  fit.cv$model %>% expect_type("character")
  fit.cv$scale %>% expect_type("logical")
  fit.cv$method %>% expect_type("character")
  fit.cv$alpha %>% expect_type("double")
  fit.cv$cv %>% expect_type("logical")
  fit.cv$fdr_control %>% expect_type("logical")
  fit.cv$call %>% expect_type("language")

  fit.cv$b %>% expect_length(dim(fit.cv$x_incidence)[2])
  fit.cv$beta %>% expect_length(dim(fit.cv$x_latency)[2])
  fit.cv$b0 %>% expect_length(1)
  fit.cv$alpha %>% expect_length(1)
  fit.cv$rate %>% expect_length(1)
  fit.cv$selected_step_inc %>% expect_length(1)
  fit.cv$selected_step_lat %>% expect_length(1)
  fit.cv$logLik %>% expect_length(1)
  fit.cv$max_c %>% expect_length(1)
  fit.cv$max.auc %>% expect_length(1)


  set.seed(4)
  temp <- generate_cure_data(n = 200, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  expect_error(cv_curegmifs(training$Time))
  training$group <- gl(2, 75)
  expect_error(cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training, x_latency = training, thresh = -1))
  expect_error(cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training, x_latency = training, epsilon = -1))
  expect_error(cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training, x_latency = training,
                         subset = group))
  expect_error(cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training, x_latency = testing))
  expect_error(cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training, x_latency = training,
                         model = "cox"))
  expect_error(cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training, x_latency = training,
                         penalty_factor_inc = penalty))
  expect_error(cv_curegmifs(Surv(Time, Censor) ~ .,
                            data = training, x_latency = training,
                            fdr_control = TRUE, fdr = 1.2))

  set.seed(16)
  temp <- generate_cure_data(n = 100, j = 15, n_true = 3, a = 1.8, rho = 0.2)
  training <- temp$training
  fit.cv <- cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training,
                         x_latency = training, fdr_control = FALSE,
                         maxit = 450, epsilon = 0.01, n_folds = 2,
                         seed = 23, verbose = TRUE
  )
  expect_equal(round(fit.cv$rate, 6), 4.518993)
  expect_equal(round(fit.cv$alpha, 7), 1.0191001)
})
