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
#'
#' @examples
#' dat <- matrix(c(1, 1, 1, 3, 0, 0, 2, 1, 0), 3, byrow = TRUE)
#' gom(dat)
gmo <- function (dat, corstr = "independence", link = "glogit", ...) {
  if (!is.matrix(dat))
    stop("dat must be a matrix")
  if (!(corstr %in% c("independence", "exchangeable")))
    stop("corstr must be one of: independence, exchangeable")
  if (!(link %in% c("glogit", "identity")))
    stop("link must be one of: glogit, identity")

  RHO <- (corstr == "exchangeable")
  ndim <- dim(dat)[[2]]
  out <- list(ndim = ndim)

  estimates <- iterative_function(dat = dat, RHO = RHO)
  out$num.iter <- estimates$num.iter
  if (RHO) out$rho <- estimates$rho

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
  # out <- iterative_function(dat = dat, RHO = RHO, ...)
  # p.vec <- unname(out$estimates[1:(ndim - 1)])
  # rho <- unname(out$estimates[ndim + 1])
  # out$var.mat <- get_var_ident(p.vec, rho, dat)
  # out$conf.int <- get_conf_int(as.matrix(p.vec[1:(ndim - 1)]),
  #                            as.matrix(sqrt(diag(out$var.mat))))
  # rownames(out$conf.int) <- c(sprintf("P%d", seq(1:(ndim - 1))))
  # colnames(out$conf.int) <- c("Lower", "Upper")
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
#'
#' @examples
#' dat <- matrix(c(1, 1, 1, 3, 0, 0, 2, 1, 0), 3, byrow = TRUE)
#' res <- gom(dat, link = "glogit")
#' print(res)
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
#'
#' @examples
#' dat <- matrix(c(1, 1, 1, 3, 0, 0, 2, 1, 0), 3, byrow = TRUE)
#' res <- gom(dat, link = "identity")
#' print(res)
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
