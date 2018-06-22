bi <- function (gmo, const.vec = NULL) {
  if (!(class(gmo) == "gmo.ident"))
    stop('gmo must be calculated using link=identity')
  m <- gmo$ndim - 1
  if (m > 1 && !(length(const.vec) == m))
    stop('constant vector must be same dimension as order of multinom')

  out <- list()

  if (m == 1) {
    out$bi <- const.vec %*% gmo$p.vec - 1
    out$var <- 4*gmo$var.mat
  } else {
    out$bi <- const.vec %*% gmo$p.vec[1:m]
    out$var <- get_bi_var(gmo$var.mat, const.vec)
  }

  out$ci <- get_conf_int(out$bi, sqrt(out$var), gmo$num.clus - 1)
  class(out) <- "bi"

  out
}

print.bi <- function (x, ...) {
  cat("BI: \n")
  print(round(x$bi, 6), ...)
  cat("\nVariance: \n")
  print(round(x$var, 6), ...)
  cat("\nConfidence Interval: \n")
  print(round(x$ci, 6), ...)
}
