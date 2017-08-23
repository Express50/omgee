library(multinom)
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

      gee[[i]] <- gom(dat, corstr = "exchangeable", ...)
      ci <- gee[[i]]$conf.int
      gee.cp[i, ] <- ci[, 1] <= p.vec[1:mord] & p.vec[1:mord] <= ci[, 2]
      gee.naive[[i]] <- gom(dat, ...)
  }

  return(list(
    gee = gee,
    gee.cp = gee.cp,
    gee.naive = gee.naive,
    gee.naive.cp = gee.naive.cp
  ))
}

res <- RunSimulation(num.sim = 500,
                     p.vec = c(0.3, 0.2, 0.4),
                     # p.vec = c(0.5, 0.5),
                     rho = 0.7,
                     num.clus = 25)
# head(res$gee)
# head(res$gee.naive)
head(res$gee.cp)
colMeans(res$gee.cp)
