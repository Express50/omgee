#' Get GEE identity Co-Variance Matrix
#'
#' Calculates GEE identity variance co-variance matrix.
#'
#' @param p.vec vector of probabilities
#' @param rho overdispersion parameter
#' @param dat matrix of multinomial clusters
#'
#' @return matrix of covariance
get_var_ident <- function (p.vec, rho, dat) {
  n <- dim(dat)[[2]]  # dimensions of multinom
  m <- n - 1  # order of multinom
  ni <- rowSums(dat) # cluster sizes
  num.clus <- dim(dat)[[1]]

  # helper vars
  rho.sq <- rho ^ 2
  dispi <- ni * (1 + (ni - 1) * rho.sq)
  p.mat <- matrix(p.vec, num.clus, m, byrow = TRUE)

  r.vec <- dat[, 1:m] - ni * p.mat
  var.mat <- diag(p.vec, m) - p.vec %*% t(p.vec)
  var.mat.inv <- solve(var.mat)

  Bi <- var.mat.inv * sum(ni ^ 2 / dispi)
  Mi <- diag(numeric(m), m)

  for (i in 1:num.clus) {
    ri <- r.vec[i, ]
    Yi.var <- ri %*% t(ri)
    Vi.inv <- (1 / dispi[i]) * var.mat.inv
    Di <- ni[i] * diag(numeric(m) + 1, m)
    Mi <- Mi + t(Di) %*% Vi.inv %*% Yi.var %*% Vi.inv %*% Di
  }

  cov.mat <- ( (sum(ni) - 1) / (sum(ni) - m) ) * (num.clus / (num.clus - 1)) * (solve(Bi) %*% Mi %*% solve(Bi))
  cov.mat
}

#' Get GEE glogit Co-Variance Matrix
#'
#' Calculates GEE glogit variance co-variance matrix.
#'
#' @param betas vector of betas
#' @param rho overdispersion parameter
#' @param dat matrix of multinomial clusters
#'
#' @return matrix of covariance
get_var_glogit <- function(betas, rho, dat) {
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
  bi <- matrix(0, m, m)
  mi <- diag(numeric(m), m)
  var.mat <- diag( p.vec, m ) - p.vec %*% t( p.vec )
  var.mat.inv <- solve(var.mat)
  r.vec <- dat[, 1:m] - ni * p.mat
  d.bar <- create_beta_equations(betas, rho, dat) / num.clus

  for (i in 1:num.clus) {
    ri <- r.vec[i, ]
    vi.inv <- (1 / dispi[i]) * var.mat.inv
    di <- ni[i] * (diag( betas.exp * denom, m ) - betas.exp %*% t( betas.exp )) / (denom ^ 2)
    dj <- t(di) %*% vi.inv %*% ri
    bi <- bi + t(di) %*% vi.inv %*% di
    mi <- mi + (dj - d.bar) %*% t(dj - d.bar)
  }

  cov.mat <- ( (sum(ni) - 1) / (sum(ni) - m) ) * (num.clus / (num.clus - 1)) * ( solve(bi) %*% mi %*% solve(bi) )
  cov.mat
}
