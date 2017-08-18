context("Iteratively solving GEEs")

dat <- matrix(
  c(18, 27, 75, 25, 184, 70, 19, 14, 18),
  nrow = 3,
  ncol = 3
)

test_that("runs with RHO=FALSE", {
  expected <- c(0.27, 0.62, 0.11)
  actual <- multinom::IterativeFunction(1, c(0.00000001, 0.00000001, 0.00000001), dat = dat, RHO=FALSE)
  expect_equal(unname(actual)[1:3], expected, tolerance = 0.01)
})

test_that("runs with RHO=TRUE", {
  expected <- c(0.29, 0.55, 0.16)
  actual <- multinom::IterativeFunction(10, c(0.00000001, 0.00000001, 0.00000001), dat = dat, RHO=TRUE)
  expect_equal(unname(actual)[1:3], expected, tolerance = 0.01)
})
