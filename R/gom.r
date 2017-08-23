#' GOM
#'
#' Creates an General Overdispersion Model for the given data
#'
#'
#' @param dat matrix of multinomial clusters
#' @param corstr one of: independence, exchangeable
#' @param method one of: GEE, DM, RCM
#' @param ... parameters to be passed on to \code{IterativeFunction}
#'
#' @return a GOM object
#' @export
#'
#' @examples
#' dat <- matrix(c(1, 1, 1, 3, 0, 0, 2, 1, 0), 3, byrow = TRUE)
#' gom(dat)
gom <- function (dat, corstr = "independence", method = "GEE", ...) {
  if (!(corstr %in% c("independence", "exchangeable")))
    stop("corstr must be one of: independence, exchangeable")
  if (!(method %in% c("GEE", "DM", "RCM")))
    stop("method must be one of: GEE, DM, RCM")

  RHO <- (corstr == "exchangeable")
  ndim <- ncol(dat)
  out <- NULL
  if (method == "GEE") {
    out <- IterativeFunction(dat = dat, RHO = RHO, ...)
    p.vec <- unname(out$estimates[1:(ndim - 1)])
    rho <- unname(out$estimates[ndim + 1])
    out$var.mat <- GetIdentityVarCov(p.vec, rho, dat)
    out$conf.int <- GetConfInt(as.matrix(p.vec[1:(ndim - 1)]),
                               as.matrix(sqrt(diag(out$var.mat))))
    rownames(out$conf.int) <- c(sprintf("P%d", seq(1:(ndim - 1))))
    colnames(out$conf.int) <- c("Lower", "Upper")
    class(out) <- "gom.gee"
  }

  out
}

#' Print GOM GEE
#'
#' Prints OM GEE object
#'
#' @param gom.gee om gee object
#'
#' @return None
#' @export
#'
#' @examples
#' dat <- matrix(c(1, 1, 1, 3, 0, 0, 2, 1, 0), 3, byrow = TRUE)
#' res <- gom(dat, method = "GEE")
#' print(res)
print.gom.gee <- function(gom.gee) {
  cat('Number of Iterations: ', gom.gee$num.iter, '\n')
  cat('\nEstimates: \n')
  print(round(gom.gee$estimates, 4))
  cat('\nConfidence Intervals: \n')
  print(round(gom.gee$conf.int, 4))

  cat('\nVar-Covar Matrix: \n')
  print(round(gom.gee$var.mat, 4))
}
registerS3method("print", "gom.gee", "print.gom.gee")
