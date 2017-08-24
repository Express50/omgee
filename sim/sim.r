library(omgee)
library(OverdispersionModelsInR)

RunSimulation <- function(num.sim, p.vec, rho, num.clus, ...) {
  ndim <- length(p.vec)
  mord <- ndim - 1
  gee <- vector("list", num.sim)
  gee.cp <- matrix(0, num.sim, mord)
  gee.naive <- vector("list", num.sim)
  gee.naive.cp <- matrix(0, num.sim, mord)

  for (i in 1:num.sim) {
      if (i %% 100 == 0) cat(sprintf("Running simulation rep %d\n", i))

      # bin <- rep(msize, num.clus)
      bin <- rnbinom(num.clus, mu = 100, size = 2) + 1
      dm <- r.dm(num.clus, p.vec, rho, bin)
      # dm <- r.rcm(num.clus, p.vec, rho, bin)
      dat <- t(dm)

        # BETAS
      # gee[[i]] <- gom2(dat, corstr = "exchangeable", ...)
      # # ci <- gee[[i]]$conf.int
      # ci <- GetConfInt(gee[[i]]$estimates[1], sqrt(diag(gee[[i]]$var.mat))[1])
      # ci.p <- exp(ci) / (1 + exp(ci))
      # gee.cp[i, ] <- ci.p[1] <= p.vec[1] & p.vec[1] <= ci.p[2]

        # P
      gee[[i]] <- gom(dat, corstr = "exchangeable", ...)
      # ci <- gee[[i]]$conf.int
      ci <- get_conf_int(gee[[i]]$estimates[1], sqrt(diag(gee[[i]]$var.mat))[1])
      # ci.p <- exp(ci) / (1 + exp(ci))
      gee.cp[i, ] <- ci[1] <= p.vec[1] & p.vec[1] <= ci[2]


      gee.naive[[i]] <- gom(dat, ...)
      ci <- get_conf_int(gee.naive[[i]]$estimates[1], sqrt(diag(gee.naive[[i]]$var.mat))[1])
      ci.p <- exp(ci) / (1 + exp(ci))
      gee.naive.cp[, ] <- ci.p[1] <= p.vec[1] & p.vec[1] <= ci.p[2]

  }

  return(list(
    gee = gee,
    gee.cp = gee.cp,
    gee.naive = gee.naive,
    gee.naive.cp = gee.naive.cp
  ))
}

res <- RunSimulation(num.sim = 1000,
                     # p.vec = c(0.3, 0.2, 0.4),
                     p.vec = c(0.9, 0.1),
                     rho = 0.7,
                     num.clus = 100)
# head(res$gee)
# head(res$gee.naive)
# head(res$gee.cp)
print(colMeans(res$gee.cp))
print(colMeans(res$gee.naive.cp))

# FOR REFERENCE
# library(msm)
# deltamethod(~exp(x1)/(1+exp(x1)+exp(x2)), c(1.0885746, 0.9947765), matrix(c(0.2792501,0.2581343,0.2581343,0.2480093), 2, 2, byrow=TRUE))

