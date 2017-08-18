context("Creating GEEs")

betas <- c(1, 1)
rho <- 0.75
dat <- matrix(1, 3, 3)

test_that("creates beta equations correctly", {
  expected <- matrix(c(-624, -624), byrow = FALSE)
  actual <- multinom::CreateBetaEquations(betas = betas,
                                          rho = rho,
                                          dat = dat)

  expect_equal(actual, expected, tolerance = 0.01)
})

test_that("creates rho equations correctly", {
  expected <- -2.98
  actual <- multinom::CreateRhoEquations(rho = rho,
                                         betas = betas,
                                         dat = dat)

  expect_equal(actual, expected, tolerance = 0.01)
})
