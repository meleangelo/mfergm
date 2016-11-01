#' Computes entropy of matrix
#' 
#' this function computes the entropy of (symmetric) matrix \code{mu}. The input should be
#' a \code{n}*\code{n} symmetric matrix, with zeros in the diagonal.
#' 
#' @param mu a symmetric matrix with zeroes in the diagonal
#' 
#' @author Angelo Mele 
#' @examples 
#' mu <- matrix(c(0,0.2,0.2,0.2,0,0.2,0.2,0.2,0), 
#'            nrow = 3, 
#'            ncol = 3
#'            )
#' entropy(mu)
#' 
#' 

entropy <- function(mu) {
  n <- dim(mu)[1]
  ent <- 0
  entropy <- .Fortran("entropy",
                n = as.integer(n),
                mu = matrix(as.double(mu), nrow = n),
                ent = as.double(ent)
                )
  return(entropy$ent)
}