loglikmf.model2new <- function(theta, addpars) {
#  loglikmf.model2 <- function(theta, n, tobs, x, ninit = 50) {
    n <- addpars$n
  tobs <- addpars$tobs
  g <- addpars$net
  nnn<-5+n
  x <- addpars$x
  ninit <- addpars$ninit #[length(addpars)]
  
  # compute potential
  potential <- as.numeric(theta[1]*tobs[1] + theta[2]*tobs[2] + 0.5*theta[3]*tobs[3])
  ## mf approximationa
  # attributes
  
  
  # parameters location
  p <- c(2,3)
  q <- 1
  
  # change stats to use
  dt <- matrix(0, ncol = 2, nrow = p[2])
  dt[,1] <- c(1,1,1)
  dt[,2] <- c(1,9,1)
  
  g <- as.sociomatrix(g)
  Ixx <- x %*% t(x) + (1-x) %*% t(1-x)
  B  <- matrix(1, nrow = n, ncol =n) %*% g
  B <- B + t(B)
  #mu <- 2*theta[1] + 2*theta[2]*Ixx + theta[2]*B/n
  mu <- theta[1] + theta[2]*Ixx + theta[2]*B/n
  mu <- exp(mu)
  mu <- mu/(1+mu)
  
  
  tol <- 0.0001
  logc <- rep(0,ninit)
  for (sim in 1:ninit) {
  #  mu <- matrix(runif(n^2,0,1), nrow = n, ncol = n)
    diag(mu) <- 0 
    psi <- logconstant_mf2(p, mu, theta, x, dt, tol) 
    
    logc[sim] <- psi$psi
  }

  psi_mf <- max(logc)

  # compute log likelihood value
  loglik <- potential - psi_mf
  #loglik <- -loglik
  return(loglik)
}