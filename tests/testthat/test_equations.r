betas <- c(1, 1)
rho <- 0
data <- matrix(1, 3, 3)

test_that("creates beta equations correctly", {
  expected <- matrix(c(-209, -209), byrow = FALSE)
  actual <- multinom::CreateBetaEquations(betas, rho, data)

  expect_equal(actual, expected, tolerance = 0.01)
})

test_that("creates rho equations correctly", {
  # TODO
  expect_equal(1, 0, label = "TODO")
})
