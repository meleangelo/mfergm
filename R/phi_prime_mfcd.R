#' second derivative of variational problem in Chatterjee-Diaconis 2013
#' 
#' this function is the second derivative of the variational problem 
#' 
#' @param p prob of link (mean field approximation)
#' @param a parameter for edges term
#' @param b parameter for 2-stars term
#' 
#' @author Angelo Mele 
#' @examples 
#' phi_prime_mfcd(p = 0.2, a = -2, b = 0.2)
#' 
#' 

#phi_prime_mfcd<-function(p,a,b) {
#    phi_prime_mfcd<-2*b*exp(2*a+2*b*p)/(( 1+ exp(2*a+2*b*p) )^2)
#    return(phi_prime_mfcd)
#}
phi_prime_mfcd<-function(p,a,b) {
  #phi_prime_mfcd<- 2*b*exp(2*a+2*b*p)/(( 1+ exp(2*a+2*b*p) )^2)
  phi_prime_mfcd<- 2*b*p*(1-p)
  return(phi_prime_mfcd)
}