#' solves Chatterjee-Diaconis 2013 variational problem for model with homogeneous players and
#' only edges and 2-stars
#' 
#' this function computes solutions to Chatterjee-Diaconis 2013 variational problem.
#' Returns a matrix with solution (row 1) and first derivative at the solution (row 2)
#' 
#' @param a parameter for edges term
#' @param b parameter for 2-stars term
#' 
#' @author Angelo Mele 
#' @examples 
#' mfcd(a = -2, b = 0.2)
#' 
#' 
mfcd <- function(a, b) {    
    solution<-uniroot.all(phi_mfcd, c(0.0000001, 0.99999999), tol = 0.0001, a = a, b = b)
    derivative<-phi_prime_mfcd(solution, a, b)
    mfcd <- rbind(solution,derivative)
    return(mfcd)
}


