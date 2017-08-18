library(multinom)
library(OverdispersionModelsInR)

RunSimulation <- function(num.sim, ndim) {
  gee.est <- matrix(0, num.sim, ndim)
  gee.se <- matrix(0, num.sim, ndim - 1)
  gee.naive.est <- matrix(0, num.sim, ndim)
  gee.naive.se <- matrix(0, num.sim, ndim - 1)

  for (i in 1:num.sim) {
    if (i %% 100 == 0) cat(sprintf("Running simulation rep %d\n", i))
    p.vec <- rep(1 / ndim, ndim)
    rho <- 0.25
    num.clus <- 100

    bin <- rnbinom(num.clus, mu = 100, size = 2) + 1
    dm <- r.dm(num.clus, p.vec, rho, bin)
    dat <- t(dm)

    gee <- IterativeFunction(10, rep(0.00000001, ndim), dat, RHO = TRUE)
    gee.naive <- IterativeFunction(10, rep(0.00000001, ndim), dat, RHO = FALSE)
    gee.est[i, ] <- gee[1:ndim + 1]
    gee.se[i, ] <- GetIdentityVariance(p.vec = gee[1:ndim - 1],
                                       rho = gee[ndim:ndim + 1],
                                       dat = dat)
    gee.naive.est[i, ] <- gee.naive[1:ndim + 1]
    gee.naive.se[i, ] <- GetIdentityVariance(p.vec = gee.naive[1:ndim - 1],
                                             rho = gee.naive[ndim:ndim + 1],
                                             dat = dat)
  }

  return(list(
    round(colMeans(gee.se), 4),           #### robust (formula)
    round(sqrt(colVars(gee.est))[1:2], 4),  #### empirical (true)
    round(colMeans(gee.naive.se), 4),
    round(sqrt(colVars(gee.naive.est))[1:2], 4)
  ))
}

RunSimulation(500, 3)
