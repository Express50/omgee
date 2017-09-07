library(omgee)
library(OverdispersionModelsInR)

run_simulation <- function(num.sim, p.vec, rho, num.clus, clus.size, ...) {
  ndim <- length(p.vec)
  mord <- ndim - 1

  dm.err <- 0
  dm.est.mat <- matrix(NA, num.sim, ndim + 1)
  dm.se.mat <- matrix(NA, num.sim, ndim - 1)
  dm.ci.mat <- matrix(NA, num.sim, 3)

  gee.ind.est <- matrix(0, num.sim, ndim)
  gee.ind.se <- matrix(0, num.sim, ndim - 1)
  gee.ind.ci <- matrix(0, num.sim, 3)

  gee.exch.err <- 0
  gee.exch.est <- matrix(NA, num.sim, ndim + 1)
  gee.exch.se <- matrix(NA, num.sim, ndim - 1)
  gee.exch.ci <- matrix(NA, num.sim, 3)

  for (i in 1:num.sim) {
      if (i %% 100 == 0) cat(sprintf("Running simulation rep %d\n", i))

      # bin <- rep(msize, num.clus)
      bin <- rnbinom(num.clus, mu = 100, size = clus.size) + 1
      dm <- r.dm(num.clus, p.vec, rho, bin)
      # dm <- r.rcm(num.clus, p.vec, rho, bin)
      dat <- t(dm)

      ind <- gmo(dat, link = "identity")
      gee.ind.est[i, ] <- ind$p.vec
      ind.se <- diag(sqrt(ind$var.mat))
      ind.ci <- get_conf_int(ind$p.vec[1], ind.se[1])
      gee.ind.se[i, ] <- ind.se
      gee.ind.ci[i, ] <- c(ind.ci,
                           ind.ci[1] <= p.vec[1] && p.vec[1] <= ind.ci[2])

      exch <- gmo(dat, corstr = "exchangeable", link="identity", ...)

      if (exch$error) {
        gee.exch.err <- gee.exch.err + 1
      } else {
        exch.se <- diag(sqrt(exch$var.mat))
        exch.ci <- get_conf_int(exch$p.vec[1], exch.se[1])

        gee.exch.est[i, ] <- c(exch$p.vec, exch$rho)
        gee.exch.se[i, ] <- exch.se
        gee.exch.ci[i, ] <- c(exch.ci,
                              exch.ci[1] <= p.vec[1] && p.vec[1] <= exch.ci[2])
      }

      dm.est <- tryCatch({
        fit.dm.mle(t(dat), rowSums(dat), var.names = c(sprintf("P%d", seq(1:ndim)), "Rho"))
      },
      error = function (e) {
        warning(e)
        return(NA)
      })

      if (!is.list(dm.est) && is.na(dm.est)) {
        dm.err <- dm.err + 1
      } else {
        dm.se <- dm.est[[1]][1:(ndim - 1), 2]
        dm.ci <- get_conf_int(dm.est[[1]][, 1][1], dm.se[1])

        dm.est.mat[i, ] <- dm.est[[1]][, 1]
        dm.se.mat[i, ] <- dm.se
        dm.ci.mat[i, ] <- c(dm.ci,
                            dm.ci[1] <= p.vec[1] && p.vec[1] <= dm.ci[2])
      }
  }

  gee.exch.est <- na.omit(gee.exch.est)
  gee.exch.se <- na.omit(gee.exch.se)
  gee.exch.ci <- na.omit(gee.exch.ci)

  dm.est.mat <- na.omit(dm.est.mat)
  dm.se.mat <- na.omit(dm.se.mat)
  dm.ci.mat <- na.omit(dm.ci.mat)

  gee.ind.empse <- sqrt(colVars(gee.ind.est)[1:mord])
  gee.exch.empse <- sqrt(colVars(gee.exch.est)[1:mord])
  dm.empse <- sqrt(colVars(dm.est.mat)[1:mord])

  return(list(
    gee.ind.est = gee.ind.est,
    gee.ind.se = gee.ind.se,
    gee.ind.empse = gee.ind.empse,
    gee.ind.ci = gee.ind.ci,
    gee.exch.err = gee.exch.err,
    gee.exch.est = gee.exch.est,
    gee.exch.se = gee.exch.se,
    gee.exch.empse = gee.exch.empse,
    gee.exch.ci = gee.exch.ci,
    dm.err = dm.err,
    dm.est = dm.est.mat,
    dm.se = dm.se.mat,
    dm.empse = dm.empse,
    dm.ci = dm.ci.mat
  ))
}

res <- run_simulation(num.sim = 500,
                      # p.vec = c(0.3, 0.2, 0.4),
                      p.vec = c(0.9, 0.1),
                      rho = 0.7,
                      num.clus = 100,
                      clus.size = 2)

cat("gee.ind: \n")
print(colMeans(res$gee.ind.est))
print(colMeans(res$gee.ind.se))
print(res$gee.ind.empse)
print(mean(res$gee.ind.ci[, 3]))
#
cat("gee.exch: \n")
print(colMeans(res$gee.exch.est))
print(colMeans(res$gee.exch.se))
print(res$gee.exch.empse)
print(mean(res$gee.exch.ci[, 3]))
cat("num errors: ", res$gee.exch.err, "\n")
#
cat("dm: \n")
print(colMeans(res$dm.est))
print(colMeans(res$dm.se))
print(res$dm.empse)
print(mean(res$dm.ci[, 3]))
cat("num errors: ", res$dm.err, "\n")

# FOR REFERENCE
# library(msm)
# deltamethod(~exp(x1)/(1+exp(x1)+exp(x2)), c(1.0885746, 0.9947765), matrix(c(0.2792501,0.2581343,0.2581343,0.2480093), 2, 2, byrow=TRUE))
