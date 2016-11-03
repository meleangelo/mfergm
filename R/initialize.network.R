#' generate initial network as random erdos-renyi graph
#' 
#' this function generates the initial network, to start
#' the simulations of the ergm, as an erdos-renyi graph 
#' 
#' @param theta vector of parameters for the ergm
#' @param n size of the network
#' @param directed logical: FALSE if undirected network (default is FALSE)
#' 
#' @author Angelo Mele 
#' @examples 
#' n <- 10
#' theta <- c(-1,3)
#' g <- initilize.network(theta, n, directed = FALSE)
#' g
#' 
#' 

initialize.network <- function(theta, n, directed = FALSE) {
  net <- network(matrix(rbinom(n^2,1,exp(theta[1])/(1+exp(theta[1]) )  ), 
                 nrow = n, ncol = n), directed = FALSE)
  return(net)
}
  