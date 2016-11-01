loglikmfR.model1 <- function(theta, n, tobs, ninit = 50) {
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
    # create initial symmetric matrix mu
    mu <- matrix(runif(n^2,0,1), nrow = n, ncol = n)
    mu[lower.tri(mu)] = t(mu)[lower.tri(mu)]
    diag(mu) <- 0 
    # use xi to avoid underflow with logs
    xi<-log(mu/(1-mu))
    diag(xi) <- 0
    kmu <- c()
    # compute weighted degrees of each player
    degrj <- rowSums(mu) 
    kmu[1]<-theta[1]*sum(1/(1+exp(-xi))) +
      0.5*theta[2]*sum( degrj^2-degrj )/n +
      0.5*sum( log(1+exp(xi)) - xi/(1+exp(-xi)) )

    eps <- 10
    tol<-0.0001
    iter<-1
    #i<-1
  
      nminus1 <- n-1
      for (i in 1:nminus1) {
        i1<-i+1
        for (j in i1:n) {
          # xi update
          #xijnew <- 2*theta[1] + (theta[2]/n)* sum(1/(1+exp(-xi[j,])) + 1/(1+exp(-xi[i,])))
          #eps <-( 2*theta[1] + (0.5*theta[2]/n)* sum(1/(1+exp(-xi[j,])) + 1/(1+exp(-xi[i,]))) )*(1/(1+exp(-xijnew))-1/(1+exp(-xi[i,j]))) - 
      

          xi[i,j] <- 2*theta[1] + (theta[2]/n)* sum(1/(1+exp(-xi[j,])) + 1/(1+exp(-xi[i,])))
          xi[j,i] <- xi[i,j]
          # compute bound 
          degrj <- rowSums(1/(1+exp(-xi))) 
          kmu[iter+1]<-theta[1]*sum(1/(1+exp(-xi))) +
            0.5*theta[2]*sum( degrj^2-degrj )/n +
            0.5*sum( log(1+exp(xi)) - xi/(1+exp(-xi)) )
          
#           kmu[iter+1]<-theta[1]*sum(1/(1+exp(-xi))) +
#             .5*theta[2]*t(1/(1+exp(-xi))) %*% g %*% (1/(1+exp(-xi)))/(n-1) +
#             sum(log(1+exp(xi)) - xi/(1+exp(-xi)))
          eps<-abs(kmu[iter+1]-kmu[iter])  # compute change in bound
          if (eps<tol) break
          iter <- iter+1
      #print(iter)
      #print(eps)
      #print(kmu[iter])
      #print(kmu[iter+1])
        }
      }
    
    
    logc[sim] <- kmu[iter]/(n^2)
       
 
    }
  
  psi_mf <- max(logc)
  
  # compute log likelihood value
  loglik <- potential - psi_mf
  #loglik <- -loglik
  return(loglik)
}




  
