test_that("generate_cure_data works", {
  set.seed(1234)
  temp <- generate_cure_data(n = 200, j = 25, n_true = 5, a = 1.8, rho = 0.3,
                             itct_mean = 0.4, same_signs = TRUE)
  temp %>% expect_type("list")
  expect_setequal(names(temp), c("training", "testing", "parameters"))
  temp %>% expect_length(3)
  temp$training %>% expect_type("list")
  temp$testing %>% expect_type("list")
  temp$parameters %>% expect_type("list")
})
