library(multinom)
library(OverdispersionModelsInR)

RunSimulation <- function(num.sim, ndim) {
  gee.est <- matrix(0, num.sim, ndim + 1)
  gee.se <- matrix(0, num.sim, ndim - 1)
  gee.naive.est <- matrix(0, num.sim, ndim + 1)
  gee.naive.se <- matrix(0, num.sim, ndim - 1)

  for (i in 1:num.sim) {
    if (i %% 100 == 0) cat(sprintf("Running simulation rep %d\n", i))
    # p.vec <- rep(1 / ndim, ndim)
    p.vec <- c(0.45, 0.35, 0.20)
    rho <- 0.7
    num.clus <- 25

    bin <- rnbinom(num.clus, mu = 100, size = 2) + 1
    dm <- r.dm(num.clus, p.vec, rho, bin)
    dat <- t(dm)

    gee <- IterativeFunction(10, rep(0.00000001, ndim), dat, RHO = TRUE)
    gee.naive <- IterativeFunction(10, rep(0.00000001, ndim), dat, RHO = FALSE)
    gee.est[i, ] <- gee[1:(ndim + 1)]
    gee.se[i, ] <- sqrt(diag(GetIdentityVarCov(p.vec = gee[1:(ndim - 1)],
                                               rho = gee[ndim + 1],
                                               dat = dat)))
    gee.naive.est[i, ] <- gee.naive[1:(ndim + 1)]
    gee.naive.se[i, ] <- sqrt(diag(GetIdentityVarCov(p.vec = gee.naive[1:(ndim - 1)],
                                                     rho = gee.naive[ndim + 1],
                                                     dat = dat)))
  }

  return(list(
    gee.est,
    gee.se,
    gee.naive.est,
    gee.naive.se
  ))
}

res <- RunSimulation(500, 3)

GEE.rho=round(colMeans(res[[1]]),4)
GEE.naive=round(colMeans(res[[3]]),4)

GEE.rho.SE=round(colMeans(res[[2]]),4)
GEE.naive.SE=round(colMeans(res[[4]]),4)

GEE.rho.SE.emp=round(sqrt(colVars(res[[1]])),4)
GEE.naive.SE.emp=round(sqrt(colVars(res[[3]])),4)


GEE.rho
GEE.rho.SE
GEE.rho.SE.emp

GEE.naive
GEE.naive.SE
GEE.naive.SE.emp
