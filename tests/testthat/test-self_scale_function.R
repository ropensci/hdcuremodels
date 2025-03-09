test_that("self_scale function works", {
  output <- hdcuremodels:::self_scale(X = matrix(1:10, ncol = 2), scale = TRUE)
  expect_equal(round(output[1, 1], 6), -1.264911)
})
