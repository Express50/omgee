library(multinom)
library(OverdispersionModelsInR)

RunSimulation <- function(num.sim, dim) {
  GEE.SE <- matrix(0, num.sim, dim - 1)
  GEE.naive.SE <- matrix(0, num.sim, dim - 1)
  for (i in 1:num.sim) {
    if (i %% 100 == 0) cat(sprintf("Running simulation rep %d\n", i))
    p.vec <- rep(1 / dim, dim)
    rho <- 0.25
    num.clus <- 100

    bin <- rnbinom(num.clus, mu = 100, size = 2) + 1
    dm <- r.dm(num.clus, p.vec, rho, bin)
    dat <- t(dm)

    GEE <- IterativeFunction(10, rep(0.00000001, dim), dat, RHO = TRUE)
    # GEE.naive <- IterativeFunction(10, rep(0.00000001, dim), dat, RHO = FALSE)
    # GEE.SE[i, ] <- GetIdentityVariance(p.vec = GEE[1:dim],
    #                               rho = GEE[dim:dim + 1],
    #                               data = data)
    # GEE.naive.SE[i, ] <- GetIdentityVariance(p.vec = GEE.naive[1:dim],
    #                                    rho = GEE.naive[dim:dim + 1],
    #                                    dat = dat)
  }

  # round(colMeans(GEE.SE),4)           #### robust (formula)
  # #round(sqrt(colVars(GEE))[1:2],4)  #### empirical (true)
  #
  # round(colMeans(GEE.SE.naive), 4)
  # #round(sqrt(colVars(GEE.naive))[1:2],4)
}

RunSimulation(500, 3)
