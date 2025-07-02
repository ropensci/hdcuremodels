#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstats {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
#' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*

test_that("cv_cureem function works correctly", {
  library(survival)
  withr::local_seed(1234)
  temp <- generate_cure_data(n = 200, j = 25, n_true = 5, a = 1.8)
  training <- temp$training
  fit.cv <- cv_cureem(Surv(Time, Censor) ~ .,
                      data = training,
                      x_latency = training, fdr_control = FALSE,
                      grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
                      n_folds = 2, seed = 23, verbose = FALSE
  )
  fit.cv %>% expect_s3_class("mixturecure")
  expect_setequal(names(fit.cv), c("b0", "b", "beta", "logLik.inc", "logLik.lat",
                                   "selected_lambda_inc", "selected_lambda_lat",
                                   "max_c", "method", "model", "penalty", "cv",
                                   "y", "x_incidence", "x_latency",
                                   "scale", "call", "fdr_control"))
  fit.cv$b0 %>% expect_type("double")
  fit.cv$b %>% expect_type("double")
  fit.cv$beta %>% expect_type("double")
  fit.cv$logLik.inc %>% expect_type("double")
  fit.cv$logLik.lat %>% expect_type("double")
  fit.cv$selected_lambda_inc %>% expect_type("double")
  fit.cv$selected_lambda_lat %>% expect_type("double")
  fit.cv$max_c %>% expect_type("double")
  fit.cv$x_incidence %>% expect_type("double")
  fit.cv$x_latency %>% expect_type("double")
  fit.cv$y %>% expect_type("double")
  fit.cv$model %>% expect_type("character")
  fit.cv$scale %>% expect_type("logical")
  fit.cv$method %>% expect_type("character")
  fit.cv$penalty %>% expect_type("character")
  fit.cv$cv %>% expect_type("logical")
  fit.cv$fdr_control %>% expect_type("logical")
  fit.cv$call %>% expect_type("language")
  fit.cv$b %>% expect_length(dim(fit.cv$x_incidence)[2])
  fit.cv$beta %>% expect_length(dim(fit.cv$x_latency)[2])
  fit.cv$b0 %>% expect_length(1)

  fit.cv.exp <- cv_cureem(Surv(Time, Censor) ~ .,
                      data = training, model = "exponential",
                      x_latency = training, fdr_control = FALSE,
                      grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
                      n_folds = 2, seed = 23, verbose = FALSE
  )
  fit.cv.exp %>% expect_s3_class("mixturecure")
  expect_setequal(names(fit.cv.exp), c("b0", "b", "beta", "rate", "logLik.inc",
                                       "logLik.lat",
                                   "selected_lambda_inc", "selected_lambda_lat",
                                   "max_c", "method", "model", "penalty", "cv",
                                   "y", "x_incidence", "x_latency",
                                   "scale", "call", "fdr_control"))
  fit.cv.exp$b0 %>% expect_type("double")
  fit.cv.exp$b %>% expect_type("double")
  fit.cv.exp$beta %>% expect_type("double")
  fit.cv.exp$rate %>% expect_type("double")
  fit.cv.exp$logLik.inc %>% expect_type("double")
  fit.cv.exp$logLik.lat %>% expect_type("double")
  fit.cv.exp$selected_lambda_inc %>% expect_type("double")
  fit.cv.exp$selected_lambda_lat %>% expect_type("double")
  fit.cv.exp$max_c %>% expect_type("double")
  fit.cv.exp$x_incidence %>% expect_type("double")
  fit.cv.exp$x_latency %>% expect_type("double")
  fit.cv.exp$y %>% expect_type("double")
  fit.cv.exp$model %>% expect_type("character")
  fit.cv.exp$scale %>% expect_type("logical")
  fit.cv.exp$method %>% expect_type("character")
  fit.cv.exp$penalty %>% expect_type("character")
  fit.cv.exp$cv %>% expect_type("logical")
  fit.cv.exp$fdr_control %>% expect_type("logical")
  fit.cv.exp$call %>% expect_type("language")
  fit.cv.exp$b %>% expect_length(dim(fit.cv$x_incidence)[2])
  fit.cv.exp$beta %>% expect_length(dim(fit.cv$x_latency)[2])
  fit.cv.exp$b0 %>% expect_length(1)
  fit.cv.exp$rate %>% expect_length(1)

  fit.cv.wei <- cv_cureem(Surv(Time, Censor) ~ .,
                          data = training, model = "weibull", measure_inc = "auc",
                          x_latency = training, fdr_control = FALSE,
                          grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
                          n_folds = 2, seed = 23, verbose = FALSE
  )
  fit.cv.wei %>% expect_s3_class("mixturecure")
  expect_setequal(names(fit.cv.wei), c("b0", "b", "beta", "rate", "alpha",
                                       "logLik.inc", "logLik.lat",
                                       "selected_lambda_inc", "selected_lambda_lat",
                                       "max_c", "max.auc", "method", "model", "penalty",
                                       "cv", "y", "x_incidence", "x_latency",
                                       "scale", "call", "fdr_control"))
  fit.cv.wei$b0 %>% expect_type("double")
  fit.cv.wei$b %>% expect_type("double")
  fit.cv.wei$beta %>% expect_type("double")
  fit.cv.wei$rate %>% expect_type("double")
  fit.cv.wei$alpha %>% expect_type("double")
  fit.cv.wei$logLik.inc %>% expect_type("double")
  fit.cv.wei$logLik.lat %>% expect_type("double")
  fit.cv.wei$selected_lambda_inc %>% expect_type("double")
  fit.cv.wei$selected_lambda_lat %>% expect_type("double")
  fit.cv.wei$max_c %>% expect_type("double")
  fit.cv.wei$max.auc %>% expect_type("double")
  fit.cv.wei$x_incidence %>% expect_type("double")
  fit.cv.wei$x_latency %>% expect_type("double")
  fit.cv.wei$y %>% expect_type("double")
  fit.cv.wei$model %>% expect_type("character")
  fit.cv.wei$scale %>% expect_type("logical")
  fit.cv.wei$method %>% expect_type("character")
  fit.cv.wei$penalty %>% expect_type("character")
  fit.cv.wei$cv %>% expect_type("logical")
  fit.cv.wei$fdr_control %>% expect_type("logical")
  fit.cv.wei$call %>% expect_type("language")
  fit.cv.wei$b %>% expect_length(dim(fit.cv$x_incidence)[2])
  fit.cv.wei$beta %>% expect_length(dim(fit.cv$x_latency)[2])
  fit.cv.wei$b0 %>% expect_length(1)
  fit.cv.wei$rate %>% expect_length(1)
  fit.cv.wei$alpha %>% expect_length(1)
  fit.cv.wei$max.auc %>% expect_length(1)

  fit.cv.wei <- cv_cureem(Surv(Time, Censor) ~ .,
                          data = training, model = "weibull",
                          x_latency = training, fdr_control = FALSE,
                          grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
                          n_folds = 2, seed = 23, verbose = FALSE
  )
  fit.cv.wei %>% expect_s3_class("mixturecure")
  expect_setequal(names(fit.cv.wei), c("b0", "b", "beta", "rate", "alpha",
                                       "logLik.inc", "logLik.lat",
                                       "selected_lambda_inc", "selected_lambda_lat",
                                       "max_c", "method", "model", "penalty", "cv",
                                       "y", "x_incidence", "x_latency",
                                       "scale", "call", "fdr_control"))
  fit.cv.wei$b0 %>% expect_type("double")
  fit.cv.wei$b %>% expect_type("double")
  fit.cv.wei$beta %>% expect_type("double")
  fit.cv.wei$rate %>% expect_type("double")
  fit.cv.wei$alpha %>% expect_type("double")
  fit.cv.wei$logLik.inc %>% expect_type("double")
  fit.cv.wei$logLik.lat %>% expect_type("double")
  fit.cv.wei$selected_lambda_inc %>% expect_type("double")
  fit.cv.wei$selected_lambda_lat %>% expect_type("double")
  fit.cv.wei$max_c %>% expect_type("double")
  fit.cv.wei$x_incidence %>% expect_type("double")
  fit.cv.wei$x_latency %>% expect_type("double")
  fit.cv.wei$y %>% expect_type("double")
  fit.cv.wei$model %>% expect_type("character")
  fit.cv.wei$scale %>% expect_type("logical")
  fit.cv.wei$method %>% expect_type("character")
  fit.cv.wei$penalty %>% expect_type("character")
  fit.cv.wei$cv %>% expect_type("logical")
  fit.cv.wei$fdr_control %>% expect_type("logical")
  fit.cv.wei$call %>% expect_type("language")
  fit.cv.wei$b %>% expect_length(dim(fit.cv$x_incidence)[2])
  fit.cv.wei$beta %>% expect_length(dim(fit.cv$x_latency)[2])
  fit.cv.wei$b0 %>% expect_length(1)
  fit.cv.wei$rate %>% expect_length(1)
  fit.cv.wei$alpha %>% expect_length(1)

  fit.cv.fdr <- cv_cureem(Surv(Time, Censor) ~ .,
                          data = training,
                          x_latency = training, model = "weibull", penalty = "lasso",
                          fdr_control = TRUE, grid_tuning = FALSE, nlambda_inc = 10,
                          nlambda_lat = 10, n_folds = 2, seed = 23, verbose = TRUE
  )
  fit.cv.fdr %>% expect_s3_class("mixturecure")
  expect_setequal(names(fit.cv.fdr), c("b0", "b", "beta", "rate", "alpha",
                                   "selected_index_inc", "selected_index_lat",
                                   "method", "model", "penalty", "cv",
                                   "y", "x_incidence", "x_latency",
                                   "scale", "call", "fdr_control"))
  fit.cv.fdr$b0 %>% expect_type("double")
  fit.cv.fdr$b %>% expect_type("double")
  fit.cv.fdr$beta %>% expect_type("double")
  fit.cv.fdr$rate %>% expect_type("double")
  fit.cv.fdr$alpha %>% expect_type("double")
  fit.cv$selected_lambda_inc %>% expect_type("double")
  fit.cv$selected_lambda_lat %>% expect_type("double")
  fit.cv$x_incidence %>% expect_type("double")
  fit.cv$x_latency %>% expect_type("double")
  fit.cv$y %>% expect_type("double")
  fit.cv$model %>% expect_type("character")
  fit.cv$scale %>% expect_type("logical")
  fit.cv$method %>% expect_type("character")
  fit.cv$penalty %>% expect_type("character")
  fit.cv$cv %>% expect_type("logical")
  fit.cv$fdr_control %>% expect_type("logical")
  fit.cv$call %>% expect_type("language")
  fit.cv.fdr$b %>% expect_length(dim(fit.cv.fdr$x_incidence)[2])
  fit.cv.fdr$beta %>% expect_length(dim(fit.cv.fdr$x_latency)[2])
  fit.cv.fdr$b0 %>% expect_length(1)
  fit.cv.fdr$alpha %>% expect_length(1)
  fit.cv.fdr$rate %>% expect_length(1)
  expect_equal(round(fit.cv.fdr$b0, 5), 0.32659)
  expect_equal(round(fit.cv.fdr$rate, 5), 1.71472)
  expect_equal(round(fit.cv.fdr$alpha, 6), 0.752908)

  fit.scad <- cv_cureem(Surv(Time, Censor) ~ .,
                       data = training,
                       x_latency = training, model = "cox", penalty = "SCAD",
                       fdr_control = FALSE, grid_tuning = FALSE, nlambda_inc = 10,
                       nlambda_lat = 10, n_folds = 2, seed = 23, verbose = TRUE
  )
  fit.cv %>% expect_s3_class("mixturecure")
  expect_setequal(names(fit.cv), c("b0", "b", "beta", "logLik.inc", "logLik.lat",
                                   "selected_lambda_inc", "selected_lambda_lat",
                                   "max_c", "method", "model", "penalty", "cv",
                                   "y", "x_incidence", "x_latency",
                                   "scale", "call", "fdr_control"))
  fit.scad$b0 %>% expect_type("double")
  fit.scad$b %>% expect_type("double")
  fit.scad$beta %>% expect_type("double")
  fit.scad$logLik.inc %>% expect_type("double")
  fit.scad$logLik.lat %>% expect_type("double")
  fit.scad$selected_lambda_inc %>% expect_type("double")
  fit.scad$selected_lambda_lat %>% expect_type("double")
  fit.scad$max_c %>% expect_type("double")
  fit.scad$x_incidence %>% expect_type("double")
  fit.scad$x_latency %>% expect_type("double")
  fit.scad$y %>% expect_type("double")
  fit.scad$model %>% expect_type("character")
  fit.scad$scale %>% expect_type("logical")
  fit.scad$method %>% expect_type("character")
  fit.scad$penalty %>% expect_type("character")
  fit.scad$cv %>% expect_type("logical")
  fit.scad$fdr_control %>% expect_type("logical")
  fit.scad$call %>% expect_type("language")
  fit.scad$b %>% expect_length(dim(fit.cv$x_incidence)[2])
  fit.scad$beta %>% expect_length(dim(fit.cv$x_latency)[2])
  fit.scad$b0 %>% expect_length(1)

  expect_error(cv_cureem(Time ~ ., data = training,
                         x_latency = training, model = "cox", penalty = "SCAD"))
  expect_error(cv_cureem(Censor ~ ., data = training,
                         x_latency = training, model = "cox", penalty = "SCAD"))

  expect_error(cv_cureem(training$Time))
  training$subset_group <- gl(2, 75)
  training$penalty <- rnorm(dim(training)[1])
  expect_error(cv_cureem(Surv(Time, Censor) ~ .,
                      data = training, x_latency = training,
                      subset = subset_group))
  expect_error(cv_cureem(Surv(Time, Censor) ~ .,
                      data = training, x_latency = testing))
  expect_error(cv_cureem(Surv(Time, Censor) ~ .,
                         data = training, x_latency = training,
                         subset = subset_group))
  expect_error(cv_cureem(Surv(Time, Censor) ~ .,
                      data = training, x_latency = training,
                      penalty_factor_inc = penalty))
  expect_error(cv_cureem(Surv(Time, Censor) ~ .,
                      data = training, x_latency = training,
                      lambda_inc_list = -1))
  expect_error(cv_cureem(Surv(Time, Censor) ~ .,
                      data = training, x_latency = training,
                      lambda_lat_list = -1))
  expect_error(cv_cureem(Surv(Time, Censor) ~ ., penalty = "MCP",
                      data = training, x_latency = training,
                      gamma_inc = -1))
  expect_error(cv_cureem(Surv(Time, Censor) ~ ., penalty = "MCP",
                      data = training, x_latency = training,
                      gamma_lat = -1))
  expect_error(cv_cureem(Surv(Time, Censor) ~ .,
                         data = training, x_latency = training,
                         fdr.control = TRUE, fdr = 1.2))
  expect_error(cv_cureem(Surv(Time, Censor) ~ .,
                         data = training, x_latency = training,
                         lambda_min_ratio_inc  = 1.3))
  expect_error(cv_cureem(Surv(Time, Censor) ~ .,
                         data = training, x_latency = training,
                         lambda_min_ratio_lat = -0.2))

  withr::local_seed(26)
  temp <- generate_cure_data(n = 200, j = 25, n_true = 5, a = 1.8)
  training <- temp$training
  fit.cv <- cv_cureem(Surv(Time, Censor) ~ .,
                      data = training,
                      x_latency = training, fdr_control = FALSE,
                      grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
                      n_folds = 2, seed = 23, verbose = FALSE
  )
  expect_equal(round(fit.cv$b0, 7), 0.3495074)
})
