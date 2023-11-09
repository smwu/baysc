#' Run MCMC to get posterior samples for WOLCA 
#'
#' @description
#' `run_MCMC_Rcpp_wolca` runs the Gibbs sampler MCMC algorithm using Rcpp to 
#' obtain posterior samples for the two-step unsupervised WOLCA model.
#'
#' @inheritParams run_MCMC_Rcpp
#' 
#' @details
#' A Gibbs sampler updates the parameters and variables in the following order:
#' \eqn{\pi}, `c_all`, \eqn{\theta}. Class assignments are permuted every 10 
#' iterations to encourage mixing, according to a random permutation sampler.
#' 
#' @return
#' Returns list `MCMC_out` containing:
#' \describe{
#'   \item{\code{pi_MCMC}}{Matrix of posterior samples for pi. (n_iter)xK}
#'   \item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xpxKxd}
#'   \item{\code{c_all_MCMC}}{Matrix of posterior samples for c_all. (n_iter)xn}
#' }
#'
#' @seealso [run_MCMC_Rcpp()] [post_process_wolca()] [get_estimates_wolca()] 
#' [fit_probit_wolca()] [wolca()] 
#' @importFrom gtools permute
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
#' cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
#' sampling_wt <- data_vars$sample_wt
#' 
#' # Obtain dimensions
#' n <- dim(x_mat)[1]        # Number of individuals
#' J <- dim(x_mat)[2]        # Number of exposure items
#' R <- max(apply(x_mat, 2,  # Number of exposure categories
#' function(x) length(unique(x))))  
#' # Obtain normalized weights
#' kappa <- sum(sampling_wt) / n   # Weights norm. constant
#' w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
#' 
#' # Set hyperparameters
#' K <- 30
#' alpha <- rep(1, K) / K
#' eta <- rep(1, R)
#' 
#' # First initialize OLCA params
#' OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)
#' 
#' # Then run MCMC sampling
#' MCMC_out <- run_MCMC_Rcpp_wolca(OLCA_params = OLCA_params, n_runs = 50, 
#' burn = 25, thin = 5, K = K, J = J, R = R, n = n, w_all = w_all, x_mat = x_mat, 
#' alpha = alpha, eta = eta)
#' # MCMC_out
#' 
run_MCMC_Rcpp_wolca <- function(OLCA_params, n_runs, burn, thin, K, J, R, n,
                                w_all, x_mat, alpha, eta) {
  # Number of MCMC iterations to store
  n_storage <- ceiling(n_runs / thin)  
  
  # Initialize variables
  pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
  theta_MCMC <- array(NA, dim = c(n_storage, J, K, R))
  c_all_MCMC <- z_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
  
  # Initialized values
  pi <- OLCA_params$pi
  theta <- OLCA_params$theta
  c_all <- OLCA_params$c_all
  
  # Update parameters and variables
  for (m in 1:n_runs) {
    update_pi(pi = pi, w_all = w_all, c_all = c_all, K = K, alpha = alpha)
    update_c_wolca(c_all = c_all, n = n, K = K, J = J, theta = theta,
                   x_mat = x_mat, pi = pi)
    update_theta(theta = theta, J = J, K = K, R = R, eta = eta,
                 w_all = w_all, c_all = c_all, x_mat = x_mat)

    #============== Store posterior values based on thinning  ==================
    if (m %% thin == 0) {
      m_thin <- m / thin
      pi_MCMC[m_thin, ] <- pi
      theta_MCMC[m_thin, , , ] <- theta
      c_all_MCMC[m_thin, ] <- c_all
    }
    
    #============== Relabel classes every 10 iterations to encourage mixing ====
    if (m %% 10 == 0) {
      new_order <- gtools::permute(1:K)      # New permuted labels
      new_c_all <- numeric(K)        # Class assignments with new labels
      for (k in 1:K) {
        new_c_all[c_all == k] <- new_order[k]
      }
      c_all <- new_c_all             # Relabel class assignments
      pi <- pi[new_order]            # Relabel class probabilities
      theta <- theta[, new_order, ]  # Relabel item category probabilities
      print(paste0("Iteration ", m, " completed!"))
    }
  }
  
  # Discard burn-in
  warmup <- ceiling(burn / thin)
  pi_MCMC <- pi_MCMC[-(1:warmup), ]
  theta_MCMC <- theta_MCMC[-(1:warmup), , , ]
  c_all_MCMC <- c_all_MCMC[-(1:warmup), ]
  
  MCMC_out <- list(pi_MCMC = pi_MCMC, theta_MCMC = theta_MCMC,
                   c_all_MCMC = c_all_MCMC)
  return(MCMC_out)
}