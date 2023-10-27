#' Initialize OLCA model
#'
#' @description
#' `init_OLCA` initializes priors and variables for an overfitted latent class
#' analysis (OLCA) given hyperparameters.
#'
#' @param K Number of classes
#' @param n Number of individuals
#' @param p Number of exposure items
#' @param d Number of exposure categories
#' @param alpha Vector of hyperparameters for pi. Kx1
#' @param eta Vector of hyperparameters for theta. dx1
#'
#' @return
#' Returns list `OLCA_params` containing:
#' \describe{
#'   \item{\code{pi}}{Vector parameter pi for class membership probabilities. Kx1}
#'   \item{\code{theta}}{Array parameter theta for item category probabilities. pxKxd}
#'   \item{\code{c_all}}{Vector of random initial class assignments. nx1}
#' }
#' 
#' @importFrom LaplacesDemon rdirichlet rcat
#' @export
#'
#' @examples
#' K <- 30; n <- 4000; p <- 30; d <- 4
#' alpha <- rep(1, K) / K
#' eta <- rep(1, d)
#' OLCA_params <- init_OLCA(K = K, n = n, p = p, d = d, alpha = alpha, eta = eta)
#' # OLCA_params
#' 
init_OLCA <- function(K, n, p, d, alpha, eta) {
  # Prior for pi
  pi <- c(LaplacesDemon::rdirichlet(n = 1, alpha = alpha))
  
  # Initialize class assignment, c, for individuals
  c_all <- LaplacesDemon::rcat(n = n, p = pi)
  
  # Prior for theta
  theta <- array(0, dim = c(p, K, d))
  for (j in 1:p) {
    for (k in 1:K) {
      theta[j, k, ] <- c(LaplacesDemon::rdirichlet(n = 1, alpha = eta))
    }
  }
  
  OLCA_params <- list(pi = pi, theta = theta, c_all = c_all)
  return(OLCA_params)
}