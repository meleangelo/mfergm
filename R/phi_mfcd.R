#' mean-field function for the Chatterjee-Disconis variational problem
#' 
#' 
#' this function is the same as Chatterjee-Diaconis first order conditions
#' 
#' @param p prob of link (mean field approximation)
#' @param a parameter for edges term
#' @param b parameter for 2-stars term
#' 
#' @author Angelo Mele 
#' @examples 
#' phi_mfcd(p = 0.2, a = -2, b = 0.2)
#' 
#' 

#phi_mfcd<-function(p,a,b) {
#    phi_mfcd<-exp(2*a+2*b*p)/( 1+ exp(2*a+2*b*p) ) - p
#    return(phi_mfcd)
#}
phi_mfcd<-function(p,a,b) {
  phi_mfcd<-exp(2*a+2*b*p)/( 1+ exp(2*a+2*b*p) ) - p
  return(phi_mfcd)
}
