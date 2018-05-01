bi <- function (gmo, const.vec) {
  if (!(class(gmo) == "gmo.ident"))
    stop('gmo must be calculated using link=identity')
  m <- gmo$ndim - 1
  # if (!(gmo$ndim == 3))
  #   stop("gmo must be for a trinomial dataset")

  out <- list()
  out$bi <- sum(const.vec * gmo$p.vec[1:m])
  out$var <- get_bi_var(gmo$var.mat, const.vec)
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
