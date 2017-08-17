#' Get GEE Identity Variance
#'
#' Calculates GEE identity variance.
#'
#' @param p.vec vector of probabilities
#' @param rho overdispersion parameter
#' @param data matrix of multinomial clusters
#'
#' @return standard error
#' @export
GetIdentityVariance <- function (p.vec, rho, data) {
  n <- dim(data)[[2]]  # dimensions of multinom
  m <- n - 1  # order of multinom
  ni <- rowSums(data) # cluster sizes
  num.clus <- dim(data)[[1]]

  # helper vars
  rho.sq <- rho ^ 2
  dispi = ni * (1 + (ni - 1) * rho.sq)
  p.mat <- matrix(p.vec, num.clus, m, byrow = TRUE)

  r.vec <- data[, 1:m] - ni * p.mat
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

  cov.mat <- (num.clus / (num.clus - 1)) * (solve(Bi) %*% Mi %*% solve(Bi))
  se <- sqrt(diag(cov.mat))

  se
}