bi <- function (gmo) {
  if (!(class(gmo) == "gmo.ident"))
    stop('gmo must be calculated using link=identity')
  if (!(gmo$ndim == 3))
    stop("gmo must be for a trinomial dataset")

  out <- list()
  out$bi <- gmo$p.vec[1] - gmo$p.vec[2]
  out$var <- get_bi_var(gmo$var.mat)
  out$ci <- get_conf_int(out$bi, sqrt(out$var), gmo$num.clus)

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
