#' GEE Beta Equations
#'
#' Creates GEE Equations for beta coefficients.
#'
#' @param betas vector of coefficients
#' @param rho overdispersion parameter
#' @param dat matrix of multinomial clusters
#'
#' @return estimating equations for betas
CreateBetaEquations <- function (betas, rho, dat) {
  n <- dim(dat)[[2]]  # dimensions of multinom
  m <- n - 1  # order of multinom
  ni <- rowSums(dat) # cluster sizes
  num.clus <- dim(dat)[[1]]

  # helper vars
  betas.exp <- exp(betas)
  rho.sq <- rho ^ 2
  denom <- 1 + sum(betas.exp)
  dispi <- (ni * (1 + (ni - 1) * rho.sq))
  p.vec <- as.vector(betas.exp / denom)
  p.mat <- matrix(p.vec, num.clus, m, byrow = TRUE)

  # create matrices
  d.mat <- diag( betas.exp * denom, m ) - betas.exp %*% t( betas.exp ) / (denom ^ 2)
  var.mat <- diag( p.vec, m ) - p.vec %*% t( p.vec )
  r.vec <- dat[, 1:m] - ni * p.mat

  # result
  d.mat %*% solve(var.mat) %*% rowSums( t((ni / dispi) * r.vec) )
}

#' GEE Rho Equations
#'
#' Creates GEE Equations for rho (overdispersion parameter).
#'
#' @param rho overdispersion parameter
#' @param betas vector of coefficients
#' @param dat matrix of multinomial clusters
#'
#' @return estimating equations for rho
CreateRhoEquations <- function (rho, betas, dat) {
  n <- dim(dat)[[2]]  # dimensions of multinom
  m <- n - 1  # order of multinom
  ni <- rowSums(dat) # cluster sizes
  num.clus <- dim(dat)[[1]]

  # helper vars
  betas.exp <- exp(betas)
  rho.sq <- rho ^ 2
  denom <- 1 + sum(betas.exp)
  dispi <- ni * (1 + (ni - 1) * rho.sq)
  p.vec <- as.vector(betas.exp / denom)
  p.mat <- matrix(p.vec, num.clus, m, byrow = TRUE)

  # create matrices
  var.mat <- diag(p.vec, m) - p.vec %*% t(p.vec)
  r.vec <- dat[, 1:m] - ni * p.mat

  res <- tr( r.vec %*% solve(var.mat) %*% t( (1 / dispi) * r.vec) ) - m * (num.clus - 1)
  res
}
