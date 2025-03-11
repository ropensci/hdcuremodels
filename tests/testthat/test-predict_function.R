test_that("predict function works correctly", {
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
  predict_train <- predict(fit)
  predict_train$p_uncured %>% expect_type("double")
  predict_train$linear_latency %>% expect_type("double")
  predict_train$latency_risk %>% expect_type("character")
  expect_setequal(names(predict_train), c("p_uncured", "linear_latency",
                                       "latency_risk"))
  expect_equal(round(predict_train$p_uncured[1], 7), 0.8598186)
  expect_equal(round(predict_train$linear_latency[1], 7), -1.212311)
  expect_equal(predict_train$latency_risk[1], "low risk")
  predict_test <- predict(fit, newdata = testing)
  predict_test$p_uncured %>% expect_type("double")
  predict_test$linear_latency %>% expect_type("double")
  predict_test$latency_risk %>% expect_type("character")
  fit.cv <- cv_cureem(Surv(Time, Censor) ~ .,
                      data = training,
                      x_latency = training, fdr_control = FALSE,
                      grid_tuning = FALSE, nlambda_inc = 10, nlambda_lat = 10,
                      n_folds = 2, seed = 23, verbose = TRUE
  )
  predict_train_cv <- predict(fit.cv)
  predict_train_cv$p_uncured %>% expect_type("double")
  predict_train_cv$linear_latency %>% expect_type("double")
  predict_train_cv$latency_risk %>% expect_type("character")

  fit.lm <- lm(Time ~ Censor, data = training)
  expect_error(predict.mixturecure(fit.lm), "Error: class of object must be mixturecure")
})
