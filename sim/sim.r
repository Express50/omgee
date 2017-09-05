library(omgee)
library(OverdispersionModelsInR)

run_simulation <- function(num.sim, p.vec, rho, num.clus, clus.size, ...) {
  ndim <- length(p.vec)
  mord <- ndim - 1
  # gee.ind.glogit <- vector("list", num.sim)
  # gee.ind.ident <- vector("list", num.sim)
  # # gee.exch.glogit <- vector("list", num.sim)
  # gee.exch.ident <- vector("list", num.sim)
  dm.est.mat <- matrix(0, num.sim, ndim + 1)
  dm.se.mat <- matrix(0, num.sim, ndim - 1)
  gee.ind.est <- matrix(0, num.sim, ndim)
  gee.ind.se <- matrix(0, num.sim, ndim - 1)
  gee.exch.est <- matrix(0, num.sim, ndim + 2)
  gee.exch.se <- matrix(0, num.sim, ndim - 1)

  for (i in 1:num.sim) {
      if (i %% 100 == 0) cat(sprintf("Running simulation rep %d\n", i))

      # bin <- rep(msize, num.clus)
      bin <- rnbinom(num.clus, mu = 100, size = clus.size) + 1
      dm <- r.dm(num.clus, p.vec, rho, bin)
      # dm <- r.rcm(num.clus, p.vec, rho, bin)
      dat <- t(dm)

      ind <- gmo(dat, link = "identity")
      gee.ind.est[i, ] <- ind$p.vec
      gee.ind.se[i, ] <- diag(sqrt(ind$var.mat))
      exch <- gmo(dat, corstr = "exchangeable", link="identity", ...)
      gee.exch.est[i, ] <- c(exch$p.vec, exch$rho, exch$error)
      gee.exch.se[i, ] <- diag(sqrt(exch$var.mat))

      dm.est <- fit.dm.mle(t(dat), rowSums(dat), var.names = c(sprintf("P%d", seq(1:ndim)), "Rho"))
      dm.est.mat[, ] <- dm.est[[1]][, 1]
      dm.se.mat[i, ] <- dm.est[[1]][1:(ndim - 1), 2]
  }

  # print(gee.exch.est)

  gee.ind.empse <- sqrt(colVars(gee.ind.est)[1:mord])

  return(list(
    gee.ind.est = gee.ind.est,
    gee.ind.se = gee.ind.se,
    gee.ind.empse = gee.ind.empse,
    gee.exch.est = gee.exch.est
  ))
}

res <- run_simulation(num.sim = 500,
                      # p.vec = c(0.3, 0.2, 0.4),
                      p.vec = c(0.9, 0.1),
                      rho = 0.7,
                      num.clus = 150,
                      clus.size = 2)

cat("gee.ind: \n")
print(colMeans(res$gee.ind.est))
print(colMeans(res$gee.ind.se))
print(res$gee.ind.empse)

cat("gee.exch: \n")
errors = which(res$gee.exch.est[, 4] == 1)
gee.exch.est.clean <- res$gee.exch.est
if (length(errors) > 0) gee.exch.est.clean <- res$gee.exch.est[-errors, ]
print(colMeans(gee.exch.est.clean))

# print(colMeans(gee.ind.est))
# print(colMeans(gee.ind.se))
# print(sqrt(colVars(gee.ind.est)))
#
# print(colMeans(res$dm.est))
# print(colMeans(res$dm.se))
# print(sqrt(colVars(res$dm.est)))



# FOR REFERENCE
# library(msm)
# deltamethod(~exp(x1)/(1+exp(x1)+exp(x2)), c(1.0885746, 0.9947765), matrix(c(0.2792501,0.2581343,0.2581343,0.2480093), 2, 2, byrow=TRUE))
