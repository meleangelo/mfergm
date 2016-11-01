loglikmf.model2wrap <- function(theta, addpars) {
  n <- addpars$n
  tobs <- addpars$tobs
  nnn<-5+n
  x <- addpars$x
  ninit <- addpars$ninit #[length(addpars)]
  fn <- loglikmf.model2(theta, n, tobs, x, ninit) 
  return(fn)  
}