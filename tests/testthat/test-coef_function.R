#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*

test_that("coef function works correctly", {
  library(survival)
  set.seed(1234)
  temp <- generate_cure_data(n = 100, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  fit <- curegmifs(Surv(Time, Censor) ~ .,
                   data = training, x_latency = training,
                   model = "weibull", thresh = 1e-4, maxit = 2000,
                   epsilon = 0.01, verbose = FALSE
  )
  output <- coef(fit)
  output$rate %>% expect_length(1)
  output$shape %>% expect_length(1)
  output$b0 %>% expect_length(1)
  output$beta_inc %>% expect_length(dim(fit$x_incidence)[2])
  output$beta_lat %>% expect_length(dim(fit$x_latency)[2])
  output$rate %>% expect_type("double")
  output$shape %>% expect_type("double")
  output$b0 %>% expect_type("double")
  output$beta_inc %>% expect_type("double")
  output$beta_lat %>% expect_type("double")
  expect_equal(round(output$rate, 6), 1.908932)
  expect_equal(round(output$shape, 6), 1.759779)
  expect_equal(round(output$b0, 7), 0.6942815)
  expect_error(coef("x"))
  expect_error(coef(1))
  expect_error(coef(fit, model_select = "aic"))
  expect_error(coef(fit, model_select = "caic"))
  expect_error(coef(fit, model_select = "maic"))
  expect_error(coef(fit, model_select = "bic"))
  expect_error(coef(fit, model_select = "mbic"))
  expect_error(coef(fit, model_select = "ebic"))
  expect_setequal(names(output), c("rate", "shape", "b0", "beta_inc", "beta_lat"))

  fit.lm <- lm(Time ~ Censor, data = training)
  expect_error(coef.mixturecure(fit.lm), "Error: class of object must be mixturecure")

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
  output <- coef(fit.cv)
  output$b0 %>% expect_length(1)
  output$beta_inc %>% expect_length(dim(fit.cv$x_incidence)[2])
  output$beta_lat %>% expect_length(dim(fit.cv$x_latency)[2])
  expect_setequal(names(output), c("b0", "beta_inc", "beta_lat"))
  output$b0 %>% expect_type("double")
  output$beta_inc %>% expect_type("double")
  output$beta_lat %>% expect_type("double")
  output$beta_inc %>% expect_vector()
  output$beta_lat %>% expect_vector()

  fit.cv.fdr <- cv_cureem(Surv(Time, Censor) ~ .,
                          data = training,
                          x_latency = training, model = "weibull", penalty = "lasso",
                          fdr_control = TRUE, grid_tuning = FALSE, nlambda_inc = 10,
                          nlambda_lat = 10, n_folds = 2, seed = 23, verbose = TRUE
  )
  output <- coef(fit.cv.fdr)
  output$rate %>% expect_length(1)
  output$shape %>% expect_length(1)
  output$b0 %>% expect_length(1)
  output$beta_inc %>% expect_length(dim(fit.cv.fdr$x_incidence)[2])
  output$beta_lat %>% expect_length(dim(fit.cv.fdr$x_latency)[2])
  expect_setequal(names(output), c("rate","shape", "b0", "beta_inc", "beta_lat"))
  output$rate %>% expect_type("double")
  output$shape %>% expect_type("double")
  output$b0 %>% expect_type("double")
  output$beta_inc %>% expect_type("double")
  output$beta_lat %>% expect_type("double")
  output$beta_inc %>% expect_vector()
  output$beta_lat %>% expect_vector()

  fit.cv.fdr <- cv_cureem(Surv(Time, Censor) ~ .,
                          data = training, x_latency = training,
                          model = "exponential", penalty = "lasso",
                          fdr_control = TRUE, grid_tuning = FALSE, nlambda_inc = 10,
                          nlambda_lat = 10, n_folds = 2, seed = 23, verbose = TRUE
  )
  output <- coef(fit.cv.fdr)
  output$rate %>% expect_length(1)
  output$b0 %>% expect_length(1)
  output$beta_inc %>% expect_length(dim(fit.cv.fdr$x_incidence)[2])
  output$beta_lat %>% expect_length(dim(fit.cv.fdr$x_latency)[2])
  expect_setequal(names(output), c("rate", "b0", "beta_inc", "beta_lat"))
  output$rate %>% expect_type("double")
  output$b0 %>% expect_type("double")
  output$beta_inc %>% expect_type("double")
  output$beta_lat %>% expect_type("double")
  output$beta_inc %>% expect_vector()
  output$beta_lat %>% expect_vector()

  fit.cv.gmifs <- cv_curegmifs(Surv(Time, Censor) ~ .,
                         data = training, model = "exponential",
                         x_latency = training, fdr_control = FALSE,
                         maxit = 450, epsilon = 0.01, n_folds = 2,
                         seed = 23, verbose = TRUE
  )
  output <- coef(fit.cv.gmifs)
  output$rate %>% expect_length(1)
  output$b0 %>% expect_length(1)
  output$beta_inc %>% expect_length(dim(fit.cv.gmifs$x_incidence)[2])
  output$beta_lat %>% expect_length(dim(fit.cv.gmifs$x_latency)[2])
  expect_setequal(names(output), c("rate", "b0", "beta_inc", "beta_lat"))
  output$rate %>% expect_type("double")
  output$b0 %>% expect_type("double")
  output$beta_inc %>% expect_type("double")
  output$beta_lat %>% expect_type("double")
  output$beta_inc %>% expect_vector()
  output$beta_lat %>% expect_vector()

  fit.cv.gmifs <- cv_curegmifs(Surv(Time, Censor) ~ .,
                               data = training, model = "weibull",
                               x_latency = training, fdr_control = TRUE,
                               maxit = 450, epsilon = 0.01, n_folds = 2,
                               seed = 23, verbose = TRUE
  )
  output <- coef(fit.cv.gmifs, model_select = "cAIC")
  output$rate %>% expect_length(1)
  output$shape %>% expect_length(1)
  output$b0 %>% expect_length(1)
  output$beta_inc %>% expect_length(dim(fit.cv.gmifs$x_incidence)[2])
  output$beta_lat %>% expect_length(dim(fit.cv.gmifs$x_latency)[2])
  expect_setequal(names(output), c("rate", "shape", "b0", "beta_inc", "beta_lat"))
  output$rate %>% expect_type("double")
  output$shape %>% expect_type("double")
  output$b0 %>% expect_type("double")
  output$beta_inc %>% expect_type("double")
  output$beta_lat %>% expect_type("double")
  output$beta_inc %>% expect_vector()
  output$beta_lat %>% expect_vector()

  fit.cv <- cv_curegmifs(Surv(Time, Censor) ~ ., data = training,
                         penalty_factor_inc = rep(c(0, 1), c(1, 11)),
                         measure_inc = "auc",
                         x_latency = training, fdr_control = FALSE,
                         maxit = 450, epsilon = 0.01, n_folds = 2,
                         seed = 23, verbose = FALSE, parallel = FALSE
  )
  output <- coef(fit.cv.gmifs, model_select = 375)
  output$rate %>% expect_length(1)
  output$shape %>% expect_length(1)
  output$b0 %>% expect_length(1)
  output$beta_inc %>% expect_length(dim(fit.cv.gmifs$x_incidence)[2])
  output$beta_lat %>% expect_length(dim(fit.cv.gmifs$x_latency)[2])
  expect_setequal(names(output), c("rate", "shape", "b0", "beta_inc", "beta_lat"))
  output$rate %>% expect_type("double")
  output$shape %>% expect_type("double")
  output$b0 %>% expect_type("double")
  output$beta_inc %>% expect_type("double")
  output$beta_lat %>% expect_type("double")
  output$beta_inc %>% expect_vector()
  output$beta_lat %>% expect_vector()
})
