#' Iterative GEE Function
#'
#' Iteratively finds coefficients and overdispersion parameter using
#' estimating equations.
#'
#' @param num.iter max number of iterations
#' @param thresholds vector of thresholds for betas and rho
#' @param data matrix of multinomial clusters
#' @param RHO if FALSE, rho will not be estimated
#'
#' @return matrix of parameters (probabilities and rho)
#' @export
IterativeFunction <- function(num.iter, thresholds, data, RHO=FALSE) {
  n <- dim(data)[[2]]
  m <- n - 1
  betas.prev <- numeric(m) + 1
  rho.prev <- 0
  betas <- betas.prev + 100 * thresholds[1:m]
  rho <- 100 * thresholds[n]

  iter <- 0

  # in this mode, assume rho to be 0
  if (RHO == FALSE) {
    # solve beta equations
    db2 <- multiroot(f = CreateBetaEquations,
                     start = betas.prev,
                     rho = rho,
                     data = data)
    # update betas
    betas <- db2$root
    rho <- 0
    iter <- db2$iter
  } else {
    # alternate between beta and rho equations until parameters are close enough
    # based on the given thresholds or until max number of iterations is reached
    while (abs( c(betas, rho) - c(betas.prev, rho.prev) ) > thresholds &&
           iter < num.iter) {
      betas.prev <- betas
      # solve beta equations
      db2 <- multiroot(f = CreateBetaEquations,
                       start = betas.prev,
                       rho = rho,
                       data = data)
      # update betas
      betas <- db2$root
      rho.prev <- rho
      # solve rho equations
      db3 <- uniroot(f = CreateRhoEquations,
                     interval = c(0, 1),
                     betas = betas,
                     data = data)
      # update rho
      rho <- db3$root
      iter <- iter + 1
    }
  }

  p.vec <- exp(betas) / (1 + sum(exp(betas)))
  params.mat <- c(p.vec, 1 - sum(p.vec), rho, iter)
  names(params.mat) <- c( sprintf("P%d", seq(1:n)), "Rho", "Num. Iterations" )

  params.mat
}
