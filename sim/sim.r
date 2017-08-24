library(omgee)
library(OverdispersionModelsInR)

run_simulation <- function(num.sim, p.vec, rho, num.clus, ...) {
  ndim <- length(p.vec)
  mord <- ndim - 1
  # gee.ind.glogit <- vector("list", num.sim)
  gee.ind.ident <- vector("list", num.sim)
  # gee.exch.glogit <- vector("list", num.sim)
  gee.exch.ident <- vector("list", num.sim)
  dm.est.mat <- matrix(0, num.sim, ndim + 1)
  dm.se.mat <- matrix(0, num.sim, ndim - 1)

  for (i in 1:num.sim) {
      if (i %% 100 == 0) cat(sprintf("Running simulation rep %d\n", i))

      # bin <- rep(msize, num.clus)
      bin <- rnbinom(num.clus, mu = 100, size = 2) + 1
      dm <- r.dm(num.clus, p.vec, rho, bin)
      # dm <- r.rcm(num.clus, p.vec, rho, bin)
      dat <- t(dm)

      gee.ind.ident[[i]] <- gmo(dat, link = "identity")
      gee.exch.ident[[i]] <- gmo(dat, corstr = "exchangeable", link="identity", ...)

      dm.est <- fit.dm.mle(t(dat), rowSums(dat), var.names = c(sprintf("P%d", seq(1:ndim)), "Rho"))
      dm.est.mat[, ] <- dm.est[[1]][, 1]
      dm.se.mat[i, ] <- dm.est[[1]][1:(ndim - 1), 2]
  }

  return(list(
    gee.ind.ident = gee.ind.ident,
    gee.exch.ident = gee.exch.ident,
    dm.est = dm.est.mat,
    dm.se = dm.se.mat
  ))
}

res <- run_simulation(num.sim = 500,
                      # p.vec = c(0.3, 0.2, 0.4),
                      p.vec = c(0.9, 0.1),
                      rho = 0.7,
                      num.clus = 200)

gee.ind.est <- do.call(rbind, do.call(rbind, res$gee.ind.ident)[, 3])
gee.ind.se <- do.call(rbind,
                      lapply(do.call(rbind, res$gee.ind.ident)[, 4],
                             function (x) diag(sqrt(x))))
print(colMeans(gee.ind.est))
print(colMeans(gee.ind.se))
print(sqrt(colVars(gee.ind.est)))

print(colMeans(res$dm.est))
print(colMeans(res$dm.se))
print(sqrt(colVars(res$dm.est)))



# FOR REFERENCE
# library(msm)
# deltamethod(~exp(x1)/(1+exp(x1)+exp(x2)), c(1.0885746, 0.9947765), matrix(c(0.2792501,0.2581343,0.2581343,0.2480093), 2, 2, byrow=TRUE))
