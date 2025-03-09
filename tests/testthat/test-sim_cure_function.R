test_that("multiplication works", {
  output <- hdcuremodels:::sim_cure(100, mu = 1, censor_mu = 3, reps = 1000,
                                    seed = 13)
  expect_equal(round(output[2], 6), 0.983344)
})
