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
To install the package you will need to install `Rtools` from WEBSITE HERE.


## How to use the package
