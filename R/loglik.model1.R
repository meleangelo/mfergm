loglik.model1 <- function(theta, tobs) {
  # compute potential
  potential <- as.numeric(theta[1]*tobs[1] + 0.5*theta[2]*tobs[2])
  # mean-field approx
  mus <- mfcd(theta[1],theta[2])  
  mus <-as.numeric(mus[1,])
  # compute kappa_mf
  kmu <- array(NA,length(mus))
  kmu <- theta[1]*mus + 0.5*theta[2]*(mus^2) - 0.5*(mus*log(mus)+(1-mus)*log(1-mus))
  psi_mf <- max(kmu)
  loglik <- potential - psi_mf
  #loglik <- -loglik
  return(loglik)
}
