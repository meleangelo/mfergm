#' Compute log constant of ERGM using mean-field approximation
#' 
#' This function computes the log-constant of the ERGM
#' using a mean-field approximation of the model
#' 
#' @param pp 2-dimensional array. the first entry 
#' indicates how many variables are in the direct utility. 
#' The second entry counts how many parameters.
#' 
#' @param mu matrix of initial values for mu
#' @param theta parameters of the ergm
#' @param x matrix with individual nodes attributes
#' @param dt matrix with indicator of which column of x and which change statistics to use
#' @param tol tolerance for the iterative process
#' 
#' @examples
#' n <- 10   # small network
#' x <- rnorm(n,0,1)  #attributes
#' mu <- matrix(runif(n^2,0,1), nrow = n, ncol = n) #initialize matrix mu
#' diag(mu) <- 0  #make sure mu has 0 diagonal

#' # Now let's get the variational mean-field aprpoximation
#' p <- c(1,2)  # 2 variables only
#' q <- 1   # x has only 1 column
#' dt <- matrix(0, ncol = p[2], nrow = 2) 
#' dt[,1] <- c(1,1)  # which column of x to consider 
#' dt[,2] <- c(1,1)  # which change stats to use

#' theta = c(-3.0,1.5)  # parameters
#' tol = 0.000001   # tolerance

#' psi <- logconstant_mf(p, mu, theta, x, dt, tol)  #compute log-constant
#' psi$psi  

logconstant_mf <- function(pp, mu, theta, x, dt, tol) {
  n <- dim(mu)[1]
  q <-max(1,dim(x)[2])
  psi <- 0
  logconstant_mf <- .Fortran("psi_mf",
                             n = as.integer(n),
                             q = as.integer(q),
                             pp = as.integer(pp),
                             mu = matrix(as.double(mu), nrow = n),
                             theta = as.double(theta),
                             x = matrix(as.double(x), nrow = n, ncol = q ),
                             dt = matrix(as.integer(dt), ncol = ncol(dt), nrow = nrow(dt)),
                             tol = as.double(tol),
                             psi = as.double(psi)  )
  #return(logconstant_mf$psi)
}