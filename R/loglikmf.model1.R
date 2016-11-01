loglikmf.model1 <- function(theta, n, tobs, ninit = 50) {
  # compute potential
  potential <- as.numeric(theta[1]*tobs[1] + 0.5*theta[2]*tobs[2])
  ## mf approximationa
  # attributes
  x <- rnorm(n,0,1)
  
  # parameters location
  p <- c(1,2)
  q <- 1
  
  # change stats to use
  dt <- matrix(0, ncol = p[2], nrow = 2)
  dt[,1] <- c(1,1)
  dt[,2] <- c(1,1)
  
  tol <- 0.0001
  logc <- rep(0,ninit)
  for (sim in 1:ninit) {
    mu <- matrix(runif(n^2,0,1), nrow = n, ncol = n)
    diag(mu) <- 0 
    psi <- logconstant_mf2(p, mu, theta , x, dt, tol) 
    
    logc[sim] <- psi$psi
  }

  psi_mf <- max(logc)

  # compute log likelihood value
  loglik <- potential - psi_mf
  #loglik <- -loglik
  return(loglik)
}