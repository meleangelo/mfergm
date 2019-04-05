# mfergm: a package for variational approximations of (some) ergm models
Angelo Mele and Lingjiong Zhu

## Intro
This package provides functions for the estimation of several specifications of ERGMs using variational mean-field approximations. 
Theoretical results are contained in 
Mele, Angelo and Lingjiong Zhu (2019), "Approximate variational estimation of a model of network formation", 


## Installation
# MAC
To install the package you need to have XCode installed, including the command line. 
Make sure that your gcc-gfortran compilers are up to date and linked. 

You will also need to install `statnet`, `rootSolve` and `optimx` from CRAN.
```{r}
install.packages("statnet")
install.packages("rootSolve")
install.packages("optimx")
```

Then you can proceed to install the package `mfergm` by

```{r}
library(devtools}
install_github("meleangelo/mfergm")
```

# Windows
To install the package you will need to install `Rtools`.


## Replication of Monte Carlo

```{r}

### TABLE n=50, sparse, nsim=100

library(mfergm)
rm(list=ls())


#### monte carlo inputs
nplayers <- 50  # number of players
true_parameters <- c(-2,1,1,1) # true parameter vector
number_of_simulations <- 100  # number of simulations
number_of_initializations <- 1 # number of times to re-start from random mu

### save estimation results in data.frame
results <-  data.frame(matrix(NA, ncol = 12, nrow = number_of_simulations))

### set seed
seed <- 1977

### simulations
results <- simulate.model4(true_parameters, n = nplayers, 
                           nsims = number_of_simulations, 
                           ninit = number_of_initializations,
                           ergm = TRUE, mfergm = TRUE, mple = TRUE, 
                           sim.seed = seed)$estimates


# make table to export in latex
table_results <- matrix(NA, ncol = 12, nrow = 5)
#table_results[1,2:13] <- colMeans(results)
for (i in 1:12) {
  table_results[1,i] <- mean(results[,i], na.rm = T)
  table_results[2,i] <- median(results[,i], na.rm = T)
  table_results[3,i] <- sd(results[,i], na.rm = T)
  table_results[4,i] <- quantile(results[,i], 0.05, na.rm=T)
  table_results[5,i] <- quantile(results[,i], 0.95, na.rm=T)
}
colnames(table_results) <- c("MCMLE", "MCMLE", "MCMLE", "MCMLE", 
                             "MF", "MF", "MF", "MF", 
                             "MPLE", "MPLE", "MPLE", "MPLE")

rownames(table_results) <- c("mean", "median", "se", "0.05", "0.95")

# latex table
### generate file name for table to be saved
tablename <- paste("table_n", nplayers, 
                   "_nsims", number_of_simulations, 
                   "_pars_", true_parameters[1], "_", 
                   true_parameters[2], "_", true_parameters[3], 
                   "_", true_parameters[4], 
                   "_2startriangles.tex",sep="" )


sink(tablename)
library(xtable)
xtable(table_results, digits = 3, 
       caption = paste("true=(", 
                       true_parameters[1],",",  
                       true_parameters[2],",",  
                       true_parameters[3],",",  
                       true_parameters[4],"); n=", nplayers, sep=""))
sink()


### generate file name for data to be saved
filename <- paste("table_n", nplayers, 
                  "_nsims", number_of_simulations, 
                  "_pars_", true_parameters[1], "_", 
                  true_parameters[2], "_", true_parameters[3],  
                  "_", true_parameters[4], 
                  "_2startriangles.RData",sep="" )

# save simulation results data for later
save.image(file = filename)


```
