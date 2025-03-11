#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*

test_that("sufficient_fu_test works correctly", {
  library(survival)
  set.seed(1234)
  temp <- generate_cure_data(n = 200, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  training$group <- gl(2, 75)
  km_fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
  output <- sufficient_fu_test(km_fit)
  expect_equal(output$n_n, 71)
  expect_equal(output$N, 150)
  output %>% expect_type("list")
  output %>% expect_length(3)
  km_fit_group <- survfit(Surv(Time, Censor) ~ group, data = training)
  output <- sufficient_fu_test(km_fit_group)
  expect_equal(output$n_n[1], 38)
  output %>% expect_type("list")
  output %>% expect_length(4)

  fit.lm <- lm(Time ~ Censor, data = training)
  expect_error(sufficient_fu_test(fit.lm), "Error: object must be a survfit object")
})
