#' GMO - General Multinomial Overdisperion Model
#'
#' Creates an General Multinomial Overdispersion Model for the given data
#'
#'
#' @param dat matrix of multinomial clusters
#' @param corstr one of: independence, exchangeable
#' @param link one of: glogit, identity
#' @param ... parameters to be passed on to \code{iterative_function}
#'
#' @return a GMO object
#' @export
gmo <- function (dat, corstr = "independence", link = "glogit", ...) {
  if (!is.matrix(dat))
    stop("dat must be a matrix")
  if (!(corstr %in% c("independence", "exchangeable")))
    stop("corstr must be one of: independence, exchangeable")
  if (!(link %in% c("glogit", "identity")))
    stop("link must be one of: glogit, identity")

  RHO <- (corstr == "exchangeable")
  ndim <- dim(dat)[[2]]
  num.clus <- dim(dat)[[1]]
  out <- list(ndim = ndim, num.clus = num.clus)

  estimates <- iterative_function(dat = dat, RHO = RHO)
  out$num.iter <- estimates$num.iter
  out$error <- estimates$error ## TEMP: FOR TESTING
  if (RHO) {
    out$rho <- estimates$rho
  }

  if (link == "glogit") {
    # return betas as estimates
    out$betas <- estimates$betas
    out$var.mat <- get_var_glogit(estimates$betas, estimates$rho, dat)

    class(out) <- "gmo.glogit"
  } else {
    # return p.vec as estimates
    p.vec <- exp(estimates$betas) / (1 + sum(exp(estimates$betas)))

    out$p.vec <- c(p.vec, 1 - sum(p.vec))
    out$var.mat <- get_var_ident(p.vec, estimates$rho, dat)

    class(out) <- "gmo.ident"
  }
  out
}

#' Print GMO glogit
#'
#' Prints GMO glogit object
#'
#' @param x gom glogit object
#' @param ... additional params to \code{print(x, ...)}
#'
#' @return None
#' @export
print.gmo.glogit <- function(x, ...) {
  cat('Number of Iterations: ', x$num.iter, '\n')

  cat('\nEstimates: \n')
  names(x$betas) <- c(sprintf("Beta%d", seq(1:length(x$betas))))
  print(round(x$betas, 6), ...)

  if (exists("rho", where = x)) {
    names(x$rho) <- "Rho"
    print(round(x$rho, 6), ...)
  }

  cat('\nVar-Covar Matrix: \n')
  print(round(x$var.mat, 6), ...)
}

#' Print GMO identity
#'
#' Prints GMO identity object
#'
#' @param x gom identity object
#' @param ... additional params to \code{print(x, ...)}
#'
#' @return None
#' @export
print.gmo.ident <- function(x, ...) {
  cat('Number of Iterations: ', x$num.iter, '\n')

  cat('\nEstimates: \n')
  names(x$p.vec) <- c(sprintf("P%d", seq(1:length(x$p.vec))))
  print(round(x$p.vec, 6), ...)

  if (exists("rho", where = x)) {
    names(x$rho) <- "Rho"
    print(round(x$rho, 6), ...)
  }

  cat('\nVar-Covar Matrix: \n')
  print(round(x$var.mat, 6), ...)
}

registerS3method("print", "gmo.glogit", "print.gmo.glogit")
registerS3method("print", "gmo.ident", "print.gmo.ident")
