#' Initialize OLCA model
#'
#' @description
#' `init_OLCA` initializes priors and variables for an overfitted latent class
#' analysis (OLCA) given hyperparameters.
#'
#' @param K Number of classes
#' @param n Number of individuals
#' @param J Number of exposure items
#' @param R Maximum number of exposure categories
#' @param alpha Kx1 vector of hyperparameters for prior for class membership 
#' probabilities \eqn{\pi}
#' @param eta JxR matrix of hyperparameters for prior for item consumption
#' level probabilities \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class 
#' \eqn{k}, assumed to be the same across classes.
#'
#' @details
#' First, class membership probabilities \eqn{\pi} are initialized by drawing 
#' from a Dirichlet distribution with hyperparameter \eqn{\alpha = 1/K} for each 
#' of the K components. Then, class assignments `c_all` are initialized for each 
#' individual by drawing from a Categorical distribution with parameter \eqn{\pi}. 
#' Finally, item consumption level probabilities \eqn{\theta} are initialized by 
#' drawing from a Dirichlet distribution with hyperparameter \eqn{\eta}, 
#' independently for each exposure item and latent class.
#' 
#' @return
#' Returns list `OLCA_params` containing:
#' \describe{
#'   \item{\code{pi}}{Vector parameter pi for class membership probabilities. Kx1}
#'   \item{\code{theta}}{Array parameter theta for item category probabilities. JxKxR}
#'   \item{\code{c_all}}{Vector of random initial class assignments. nx1}
#' }
#' 
#' @seealso [init_probit()] [swolca()] [solca()] [wolca()]
#' 
#' @importFrom LaplacesDemon rdirichlet rcat
#' @export
#'
#' @examples
#' K <- 30; n <- 4000; J <- 30; R <- 4
#' alpha <- rep(1, K) / K
#' eta <- matrix(1, nrow = J, ncol = R)
#' OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)
#' # OLCA_params
#' 
init_OLCA <- function(K, n, J, R, alpha, eta) {
  # Prior for pi
  pi <- c(LaplacesDemon::rdirichlet(n = 1, alpha = alpha))
  
  # Initialize class assignment, c, for individuals
  c_all <- LaplacesDemon::rcat(n = n, p = pi)
  
  # Prior for theta
  theta <- array(0, dim = c(J, K, R))
  for (j in 1:J) {
    for (k in 1:K) {
      theta[j, k, ] <- c(LaplacesDemon::rdirichlet(n = 1, alpha = eta[j, ]))
    }
  }
  
  OLCA_params <- list(pi = pi, theta = theta, c_all = c_all)
  return(OLCA_params)
}