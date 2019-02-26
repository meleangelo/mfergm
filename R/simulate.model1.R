#' Generate networks and estimates parameters 
#' 
#' This function simulates \code{nsims} networks, and 
#' compares estimates using 
#' the ergm routines, the Chatterje-Diaconis 2013 mean-field
#' and the more general mean-field approximation. The model
#' is a model with edges and 2-stars.
#' 
#' @param theta true parameters vector
#' @param n size of the network
#' @param nsims number of simulated networks to estimate
#' @param ninit number of trials for mean field approximation (if mfergm = TRUE)
#' @param ergm logical, if FALSE will not perform ergm estimation
#' @param cd logical, if FALSE will not perform Chatterjee-Diaconis estimation
#' @param mfergm logical, if FALSE will not perform mean-field estimation
#' @param mple logical, if FALSE will not perform MPLE estimation
#' @param sim.seed seed for the random number generator

simulate.model1 <- function(theta, n = 10, nsims = 2, ninit = 5,
                            ergm = FALSE, cd = TRUE, mfergm = FALSE, mple = FALSE,
                            sim.seed = 1977) {
  
  library(mfergm)
  set.seed(sim.seed) # seed for the simulations
  #n <- 100  # size of the network
  #nsims <- 2 # number of simulations
  x <- rnorm(n,0,1) # attributes
  p <- c(1,2) # parameters location
  q <- 1  # columns of x
  
  # change stats to use
  dt <- matrix(0, ncol = p[2], nrow = 2)
  dt[,1] <- c(1,1)
  dt[,2] <- c(1,1)
  
  #theta = c(-2.0,3.0) # true parameters
  tol = 0.000001 # tolerance for mf
  
  # create network from random erdos renyi
  cat("initializing network for simulations")
  g <- initialize.network(theta, n, directed = FALSE)
  
  # simulate nsims networks
  cat("\n")
  cat(paste("generate ", nsims, " networks", "\n"))
  g0 <- simulate(~edges+kstar(2), 
                 nsim = nsims, 
                 coef = theta*c(2,1/n ),
                 basis = g,
                 control=control.simulate(
                   MCMC.burnin=10000000,
                   MCMC.interval=10000)
  )
  
  # observed statistics for the sample
  tobs <- data.frame(matrix(NA, ncol = 2, nrow = nsims))
  names(tobs) <- c("edges", "kstar(2)")
  for (i in 1:nsims) {
    formula <- g0[[i]] ~ edges + kstar(2)
    #tobs[i,] <- summary(formula)/(c((n^2)/2, n^3))  
    tobs[i,] <- summary(formula)/(c((n^2)/2, (n^3)/2 ))  
  }
  cat(paste("observed stats", tobs, "\n"))
  # initialize data.frame with estimation results
  estim.table<- data.frame(matrix(NA, nrow = nsims, ncol = 8))
  names(estim.table) <- c("ergm", "ergm", 
                          "CD2013", "CD2013", 
                          "MF", "MF", "mple", "mple")
  
  # estimate using ergm routines
  if (ergm == TRUE) {
    cat(paste("estimating using ergm package", "\n"))
    for (i in 1:nsims) {
      cat("***********************************\n")
      cat(paste("estimating sample" ,i, "\n"))
      cat("***********************************\n")
      #set.seed(1977)
      formula <- g0[[i]] ~ edges + kstar(2)
      m1ergm <- ergm(formula, estimate = "MLE",
                     control=control.ergm(
                       main.method = "Stochastic-Approximation",
                       MCMC.burnin=100000,
                       MCMC.interval=1000,
                       init = theta*c(2,1/n ))
      )
      est.params <- m1ergm$coef
      estim.table[i,1:2] <- c(est.params[1:2,1])*c(.5,n)
    }
  }

  # estimate using chatterjee-diaconis
  if (cd == TRUE) {
    cat(paste("estimating using Chatterjee-Diaconis approximation", "\n"))
    library(optimx)
    for (i in 1:nsims) {
      cat("***********************************\n")
      cat(paste("estimating sample" ,i, "\n"))
      cat("***********************************\n")
      pars <- theta
      cd.est <- optimx(pars, loglik.model1, 
                       method = "BFGS", 
                       control = list(fnscale = -1), 
                       tobs = tobs[i,])
      estim.table[i,3:4] <- c(cd.est[1:2])
      
    }
  }
  
  # estimate with mfergm
  if (mfergm == TRUE) {
    cat(paste("estimating Mean field approximation", "\n"))
    library(optimx)
    for (i in 1:nsims) {
      cat("***********************************\n")
      cat(paste("estimating sample" ,i, "\n"))
      cat("***********************************\n")
      pars <- theta
      cd.est <- optimx(pars, loglikmf.model1, 
                       method = "BFGS", 
                       control = list(fnscale = -1), 
                       n = n,
                       tobs = tobs[i,], 
                       ninit = ninit)
      estim.table[i,5:6] <- c(cd.est[1:2])
      
    }
    
  }
  
  # estimate using ergm routines
  if (mple == TRUE) {
    cat(paste("estimating MPLE using ergm package", "\n"))
    for (i in 1:nsims) {
      cat("***********************************\n")
      cat(paste("estimating sample" ,i, "\n"))
      cat("***********************************\n")
      #set.seed(1977)
      formula <- g0[[i]] ~ edges + kstar(2)
      m1ergm <- ergm(formula, estimate = "MPLE",
                     control=control.ergm(
                       MCMC.burnin=100000,
                       MCMC.interval=1000,
                       init = theta*c(2,1/n ))
      )
      est.params <- m1ergm$coef
      estim.table[i,7:8] <- c(est.params[1:2,1])*c(.5,n)
    }
  }
  
  data <- list(g0,theta,n,nsims,tol,dt, tobs, estim.table)
  names(data) <- c("sample", "params", "n", "nsims", "tol", "changestats", "stats0", "estimates")
  return(data)
  
}