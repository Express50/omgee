context("Iteratively solving GEEs")

dat <- matrix(1, 3, 3)

test_that("runs with RHO=FALSE", {
  expected <- rep(1/3, 3)
  actual <- multinom::IterativeFunction(1, c(0.0001, 0.0001, 0.0001), dat = dat, RHO=FALSE)
  expect_equal(unname(actual)[1:3], expected, tolerance = 0.01)
})

test_that("runs with RHO=TRUE", {
  expected <- rep(1/3, 3)
  actual <- multinom::IterativeFunction(10, c(0.0001, 0.0001, 0.0001), dat = dat, RHO=TRUE)
  expect_equal(unname(actual)[1:3], expected, tolerance = 0.01)
})
