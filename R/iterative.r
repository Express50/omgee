library(rootSolve)

#' Iterative GEE Function
#'
#' Iteratively finds coefficients and overdispersion parameter using
#' estimating equations.
#'
#' @param dat matrix of multinomial clusters
#' @param RHO if FALSE, rho will not be estimated
#' @param max.iter max number of iterations
#' @param thresholds vector of thresholds for betas and rho
#'
#' @return list containing estimates (probabilities and rho) and number of iterations
iterative_function <- function(dat, RHO = FALSE, max.iter = 10, thresholds = NULL) {
  n <- dim(dat)[[2]]
  if (is.null(thresholds)) thresholds <- rep(0.00000001, n)

  m <- n - 1
  betas.prev <- numeric(m) + 1
  rho.prev <- 0
  betas <- betas.prev + 100 * thresholds[1:m]
  rho <- 100 * thresholds[n]

  iter <- 0
  e <- FALSE

  # in this mode, assume rho to be 0
  if (RHO == FALSE) {
    # solve beta equations
    db2 <- multiroot(f = create_beta_equations,
                     start = betas.prev,
                     rho = rho,
                     dat = dat)
    # update betas
    betas <- db2$root
    rho <- 0
    iter <- db2$iter
  } else {
    # alternate between beta and rho equations until parameters are close enough
    # based on the given thresholds or until max number of iterations is reached
    while (abs( c(betas, rho) - c(betas.prev, rho.prev) ) > thresholds &&
           iter < max.iter) {
      betas.prev <- betas
      # solve beta equations
      db2 <- multiroot(f = create_beta_equations,
                       start = betas.prev,
                       rho = rho,
                       dat = dat)
      # update betas
      betas <- db2$root
      rho.prev <- rho
      # solve rho equations
      db3 <- tryCatch(uniroot(f = create_rho_equations,
                              interval = c(0, 1),
                              betas = betas,
                              dat = dat),
                      error = function (err) {
                        warning("rho couldn't be estimated, try using gmo(dat, corstr = \"independence\")",
                                call. = FALSE,
                                immediate. = TRUE)
                        return(NA)
                      })

      # print(db3)
      if (!is.list(db3) && is.na(db3)) {e <- TRUE; break}

      # update rho
      rho <- db3$root
      iter <- iter + 1
    }
  }

  out <- list(betas = betas, rho = rho, num.iter = iter, error = e)
  out
}
