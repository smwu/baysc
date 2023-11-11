#' Get posterior estimates for WOLCA
#'
#' @description
#' `get_estimates_wolca` obtains posterior parameter samples and estimates prior 
#' for the unsupervised WOLCA model
#'
#' @inheritParams run_MCMC_Rcpp_wolca
#' @param MCMC_out Output from `run_MCMC_Rcpp_wolca` containing `pi_MCMC`, 
#' `theta_MCMC`, and `c_all_MCMC`
#' @param post_MCMC_out output from `post_process_wolca` containing `K_med`, `pi`, 
#' and `theta`
#' 
#' @details
#' First, duplicate classes that have the same modal exposure categories
#' for all items are combined to obtain the number of unique classes, `K_red`. 
#' Parameters are then renormalized for the unique classes and posterior median 
#' estimates are computed across MCMC iterations. Using these median estimates,
#' class assignments `c_all` are derived. 
#' 
#' @return
#' Returns list `estimates` containing:
#' \describe{
#'   \item{\code{K_red}}{Number of unique classes}
#'   \item{\code{pi_red}}{Matrix of final posterior samples for pi. Mx(K_red)}
#'   \item{\code{theta_red}}{Array of final posterior samples for theta. Mxpx(K_red)xd}
#'   \item{\code{pi_med}}{Vector of posterior median estimates for pi. (K_red)x1}
#'   \item{\code{theta_med}}{Array of posterior median estimates for theta. px(K_red)xd}
#'   \item{\code{c_all}}{Vector of final individual class assignments. nx1}
#'   \item{\code{pred_class_probs}}{Matrix of individual posterior class probabilities. nx(K_red)}
#' }
#'
#' @seealso [get_estimates()] [run_MCMC_Rcpp_wolca()] [post_process_wolca()] 
#' [fit_probit_wolca()] [wolca()] 
#' @importFrom plyr aaply
#' @importFrom matrixStats logSumExp
#' @importFrom LaplacesDemon rcat
#' @importFrom stats median
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
#' cluster_id <- data_vars$cluster_id  # Cluster indicators, nx1
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
#' # Set hyperparameters for fixed sampler
#' K <- 3
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
#' 
#' # Then run post-process relabeling
#' post_MCMC_out <- post_process_wolca(MCMC_out = MCMC_out, J = J, R = R,
#' class_cutoff = 0.05)
#'
#' # Then obtain posterior estimates
#' estimates <- get_estimates_wolca(MCMC_out = MCMC_out, 
#' post_MCMC_out = post_MCMC_out, n = n, J = J, x_mat = x_mat)
#' 
get_estimates_wolca <- function(MCMC_out, post_MCMC_out, n, J, x_mat) {
  
  #============== Identify unique classes using modal exposure categories ======
  # Posterior median estimate for theta across iterations
  theta_med_temp <- apply(post_MCMC_out$theta, c(2, 3, 4), stats::median)
  # Posterior modal exposure categories for each exposure item and reduced class
  theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
  # Identify unique classes
  unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
  # Number of unique classes
  K_red <- length(unique_classes)
  
  #============== Use new classes to adjust and re-normalize posterior samples =
  # Combine duplicated classes and re-normalize pi to sum to 1
  M <- dim(post_MCMC_out$pi)[1]                # Number of iterations
  pi_red <- post_MCMC_out$pi[, unique_classes, drop = FALSE] # Initialize pi for unique classes
  if (K_red < dim(post_MCMC_out$pi)[2]) {  # Check if there are duplicated classes
    for (k in 1:K_red) {
      # Find duplicated modal theta patterns
      dupes_k <- apply(theta_modes, 2, function(x)
        identical(x,theta_modes[, unique_classes[k]]))
      # Combine class proportions for all duplicated patterns together
      pi_red[, k] <- apply(as.matrix(post_MCMC_out$pi[, dupes_k]), 1, sum)
    }
  }
  # Re-normalize to ensure pi sums to 1 for each iteration
  pi_red = pi_red / rowSums(pi_red)
  
  # Get posterior parameter samples for unique classes for theta and xi
  theta_red <- post_MCMC_out$theta[, , unique_classes, , drop = FALSE]
  theta_red <- plyr::aaply(theta_red, c(1, 2, 3), function(x) x / sum(x), 
                           .drop = FALSE) # Re-normalize
  
  #============== Posterior median estimates ===================================
  pi_med <- apply(pi_red, 2, stats::median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), stats::median, na.rm = TRUE)
  theta_med <- plyr::aaply(theta_med, c(1, 2), function(x) x / sum(x),
                           .drop = FALSE)  # Re-normalize
  
  #============== Update c using unique classes and posterior estimates ========
  c_all <- MCMC_out$c_all_MCMC[M, ]                      # Final class assignments
  pred_class_probs <- matrix(NA, nrow = n, ncol = K_red) # Posterior class membership probabilities
  log_cond_c <- matrix(NA, nrow = n, ncol = K_red)       # Individual log-likelihood for each class
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    for (k in 1:K_red) {
      # Calculate theta component of individual log-likelihood assuming class k
      log_theta_comp_k <- 0
      for (j in 1:J) {
        log_theta_comp_k <- log_theta_comp_k + log(theta_med[j, k, x_mat[i, j]])
      }
      # Individual log-likelihood for class k
      log_cond_c[i, k] <- log(pi_med[k]) + log_theta_comp_k
    }
    # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs[i, ] <- exp(log_cond_c[i, ] - matrixStats::logSumExp(log_cond_c[i, ]))
    # Update class assignment using the posterior probabilities
    c_all[i] <- LaplacesDemon::rcat(n = 1, p = pred_class_probs[i, ])
  }
  
  estimates <- list(K_red = K_red, pi_red = pi_red, theta_red = theta_red,
                   pi_med = pi_med, theta_med = theta_med, c_all = c_all,
                   pred_class_probs = pred_class_probs)
  
  return(estimates)
}