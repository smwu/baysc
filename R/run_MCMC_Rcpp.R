#' Run MCMC to get posterior samples
#'
#' @description
#' `run_MCMC_Rcpp` runs the Gibbs sampler MCMC algorithm using Rcpp to obtain 
#' posterior samples.
#'
#' @inheritParams init_probit
#' @param OLCA_params Output list from `init_OLCA()` containing `pi`, `c_all`, 
#' and `theta`
#' @param probit_params Output list from `init_probit()` containing `xi` and `z_all`
#' @param alpha Vector of hyperparameters for pi. Kx1
#' @param eta Vector of hyperparameters for theta. dx1
#' @param x_mat Categorical exposure matrix. nxp
#' @param w_all Weights normalized to sum to n. nx1
#' @param n_runs Number of MCMC iterations
#' @param burn Burn-in period
#' @param thin Thinning factor
#' @return
#' Returns list `MCMC_out` containing:
#' \describe{
#'   \item{\code{pi_MCMC}}{Matrix of posterior samples for pi. (n_iter)xK}
#'   \item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xpxKxd}
#'   \item{\code{xi_MCMC}}{Array of posterior samples for xi. (n_iter)xKxq}
#'   \item{\code{c_all_MCMC}}{Matrix of posterior samples for c_all. (n_iter)xn}
#'   \item{\code{z_all_MCMC}}{Matrix of posterior samples for z_all. (n_iter)xn}
#'   \item{\code{loglik_MCMC}}{Vector of posterior samples for log-likelihood. (n_iter)x1}
#' }
#'
#' @export
#'
#' # examples
#' 
run_MCMC_Rcpp <- function(OLCA_params, probit_params, n_runs, burn, thin, K, p, d, n, 
                          q, w_all, x_mat, y_all, V, alpha, eta, mu0, Sig0) {
  # Number of MCMC iterations to store
  n_storage <- ceiling(n_runs / thin) 
  
  # Initialize variables
  pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
  theta_MCMC <- array(NA, dim = c(n_storage, p, K, d))
  xi_MCMC <- array(NA, dim = c(n_storage, K, q))
  c_all_MCMC <- z_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
  loglik_MCMC <- numeric(n_storage)
  loglik <- numeric(n)     # Individual log-likelihood
  lin_pred <- numeric(n)   # Linear predictor, V*xi
  
  # Initialized values
  pi <- OLCA_params$pi
  theta <- OLCA_params$theta
  c_all <- as.double(OLCA_params$c_all)  # allows updating by reference in rcpparmadillo
  xi <- probit_params$xi
  z_all <- probit_params$z_all
  
  # Update parameters and variables
  for (m in 1:n_runs) {
    update_pi(pi = pi, w_all = w_all, c_all = c_all, K = K, alpha = alpha)
    update_c(c_all = c_all, n = n, K = K, p = p, theta = theta, 
             x_mat = x_mat, pi = pi, z_all = z_all, V = V, xi = xi, 
             y_all = y_all)
    update_theta(theta = theta, p = p, K = K, d = d, eta = eta, 
                 w_all = w_all, c_all = c_all, x_mat = x_mat)
    xi <- update_xi(xi = xi, n = n, K = K, w_all = w_all, c_all = c_all,
                    z_all = z_all, V = V, mu0 = mu0, Sig0 = Sig0)
    z_all <- update_z(z_all = z_all, n = n, V = V, xi = xi, c_all = c_all,
                      y_all = y_all)
    update_loglik(loglik = loglik, n = n, p = p, c_all = c_all, 
                  theta = theta, x_mat = x_mat, pi = pi, 
                  z_all = z_all, V = V, xi = xi, y_all = y_all)
    
    #============== Store posterior values based on thinning  ==================
    if (m %% thin == 0) {
      m_thin <- m / thin
      pi_MCMC[m_thin, ] <- pi
      theta_MCMC[m_thin, , , ] <- theta
      xi_MCMC[m_thin, , ] <- xi
      c_all_MCMC[m_thin, ] <- c_all
      z_all_MCMC[m_thin, ] <- z_all
      loglik_MCMC[m_thin] <- sum(loglik)
    }
    
    #============== Relabel classes every 10 iterations to encourage mixing ====
    if (m %% 10 == 0) {
      new_order <- permute(1:K)      # New permuted labels
      new_c_all <- numeric(K)        # Class assignments with new labels
      for (k in 1:K) {
        new_c_all[c_all == k] <- new_order[k]
      }
      c_all <- new_c_all                # Relabel class assignments
      pi <- pi[new_order]               # Relabel class probabilities
      theta <- theta[, new_order, , drop = FALSE]     # Relabel item category probabilities
      xi <- xi[new_order, , drop = FALSE]  # Relabel probit coefficients
      print(paste0("Iteration ", m, " completed!"))
    }
  }
  
  # Discard burn-in
  warmup <- ceiling(burn / thin)
  pi_MCMC <- pi_MCMC[-(1:warmup), ]
  theta_MCMC <- theta_MCMC[-(1:warmup), , , ]
  xi_MCMC <- xi_MCMC[-(1:warmup), , , drop = FALSE]
  c_all_MCMC <- c_all_MCMC[-(1:warmup), ]
  z_all_MCMC <- z_all_MCMC[-(1:warmup), ]
  loglik_MCMC <- loglik_MCMC[-(1:warmup)]
  
  # Return output
  MCMC_out <- list(pi_MCMC = pi_MCMC, theta_MCMC = theta_MCMC, xi_MCMC = xi_MCMC,
                   c_all_MCMC = c_all_MCMC, z_all_MCMC = z_all_MCMC, 
                   loglik_MCMC = loglik_MCMC)
  return(MCMC_out)
}