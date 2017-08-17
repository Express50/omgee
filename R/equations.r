#' GEE Beta Equations
#'
#' Creates GEE Equations for beta coefficients.
#'
#' @param betas vector of coefficients
#' @param rho overdispersion parameter
#' @param data matrix of multinomial clusters
#'
#' @return estimating equations for betas
CreateBetaEquations <- function (betas, rho, data) {
  n <- dim(data)[[2]]  # dimensions of multinom
  m <- n - 1  # order of multinom
  ni <- rowSums(data) # cluster sizes
  num.clus <- dim(data)[[1]]

  # helper vars
  betas.exp <- exp(betas)
  rho.sq <- rho ^ 2
  denom <- 1 + sum(betas.exp)
  kC <- ni / (ni * (1 + (ni - 1) * rho.sq))
  p.vec <- as.vector(betas.exp / denom)
  p.mat <- matrix(p.vec, num.clus, m, byrow = TRUE)

  # create matrices
  d.mat <- diag( betas.exp * denom, m ) - betas.exp %*% t( betas.exp ) / (denom ^ 2)
  var.mat <- diag( p.vec, m ) - p.vec %*% t( p.vec )
  r.vec <- kC * data[, 1:m] - ni * p.mat

  # result
  d.mat %*% solve(var.mat) %*% rowSums( t(r.vec) )
}

#' GEE Rho Equations
#'
#' Creates GEE Equations for rho (overdispersion parameter).
#'
#' @param rho overdispersion parameter
#' @param betas vector of coefficients
#' @param data matrix of multinomial clusters
#'
#' @return estimating equations for rho
CreateRhoEquations <- function (rho, betas, data) {
  n <- dim(data)[[2]]  # dimensions of multinom
  m <- n - 1  # order of multinom
  ni <- rowSums(data) # cluster sizes
  num.clus <- dim(data)[[1]]

  # helper vars
  betas.exp <- exp(betas)
  rho.sq <- rho ^ 2
  denom <- 1 + sum(betas.exp)
  kC <- 1 / (ni * (1 + (ni - 1) * rho.sq))
  p.vec <- as.vector(betas.exp / denom)
  p.mat <- matrix(p.vec, num.clus, m, byrow = TRUE)

  # create matrices
  r.vec <- data[, 1:m] - ni * p.mat
  var.mat = diag(p.vec, m) - p.vec %*% t(p.vec)

  # result
  tr( r.vec %*% solve(var.mat) %*% t(kC * r.vec) ) - m * (num.clus - 1)
}
