#' Compute trace of a square matrix, that is, the sum
#' of its diagonal elements.
#'
#' @param m a square matrix
#'
#' @return trace of input matrix
tr <- function(m) {
    sum(diag(m))
}