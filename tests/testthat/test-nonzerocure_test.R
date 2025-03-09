test_that("nonzerocure_test works", {
  library(survival)
  set.seed(1234)
  temp <- generate_cure_data(n = 200, j = 10, n_true = 10, a = 1.8)
  training <- temp$training
  training$group <- gl(2, 75)
  km_fit <- survfit(Surv(Time, Censor) ~ 1, data = training)
  output <- nonzerocure_test(km_fit)
  output %>% expect_type("list")
  output %>% expect_length(4)
  expect_equal(round(output$p_value, 3), 0.039)

#  km_fit_group <- survfit(Surv(Time, Censor) ~ group, data = training)
#  output <- nonzerocure_test(km_fit_group)
#  output %>% expect_type("list")
#  output %>% expect_length(4)

})
