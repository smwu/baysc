#' Get posterior estimates
#'
#' @description
#' `get_estimates` obtains posterior parameter samples and estimates prior to
#' variance adjustment
#'
#' @inheritParams run_MCMC_Rcpp
#' @param MCMC_out Output from `run_MCMC_Rcpp` containing `pi_MCMC`, 
#' `theta_MCMC`, `xi_MCMC`, `c_all_MCMC`, `z_all_MCMC`, and `loglik_MCMC`
#' @param post_MCMC_out output from `post_process` containing `K_med`, `pi`, 
#' `theta`, `xi`
#' @return
#' Returns list `estimates` containing:
#' \describe{
#'   \item{\code{K_red}}{Number of unique classes}
#'   \item{\code{pi_red}}{Matrix of final posterior samples for pi. Mx(K_red)}
#'   \item{\code{theta_red}}{Array of final posterior samples for theta. Mxpx(K_red)xd}
#'   \item{\code{xi_red}}{Array of final posterior samples for xi. Mx(K_red)xq}
#'   \item{\code{pi_med}}{Vector of posterior median estimates for pi. (K_red)x1}
#'   \item{\code{theta_med}}{Array of posterior median estimates for theta. px(K_red)xd}
#'   \item{\code{xi_med}}{Matrix of posterior median estimates for xi. (K_red)xq}
#'   \item{\code{Phi_med}}{Vector of final individual outcome probabilities. nx1}
#'   \item{\code{c_all}}{Vector of final individual class assignments. nx1}
#'   \item{\code{pred_class_probs}}{Matrix of individual posterior class probabilities. nx(K_red)}
#'   \item{\code{loglik_med}}{Vector of final indiviudal log-likehoods. nx1} 
#' }
#'
#' @importFrom plyr aaply
#' @importFrom matrixStats logSumExp
#' @importFrom LaplacesDemon rcat
#' @importFrom stats dnorm median pnorm
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' cluster_id <- data_vars$cluster_id  # Cluster indicators, nx1
#' sampling_wt <- data_vars$sample_wt
#' 
#' # Obtain dimensions
#' n <- dim(x_mat)[1]        # Number of individuals
#' p <- dim(x_mat)[2]        # Number of exposure items
#' d <- max(apply(x_mat, 2,  # Number of exposure categories
#' function(x) length(unique(x))))  
#' # Obtain normalized weights
#' kappa <- sum(sampling_wt) / n   # Weights norm. constant
#' w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
#' 
#' # Probit model only includes latent class
#' V <- matrix(1, nrow = n)  
#' q <- ncol(V)   # Number of regression covariates excluding class assignment
#' 
#' # Set hyperparameters for fixed sampler
#' K <- 3
#' alpha <- rep(1, K) / K
#' eta <- rep(1, d)
#' mu0 <- Sig0 <- vector("list", K)
#' for (k in 1:K) {
#'   # MVN(0,1) hyperprior for prior mean of xi
#'   mu0[[k]] <- stats::rnorm(n = q)
#'   # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
#'   # components and mean variance 2.5 for a weakly informative prior on xi
#'   Sig0[[k]] <- diag(LaplacesDemon::rinvgamma(n = q, shape = 3.5, scale = 6.25), 
#'   nrow = q, ncol = q)
#' }
#' 
#' # First initialize OLCA params
#' OLCA_params <- init_OLCA(K = K, n = n, p = p, d = d, alpha = alpha, eta = eta)
#' 
#' # Then initialize probit params 
#' probit_params <- init_probit(K = K, n = n, q = q, V = V, mu0 = mu0, 
#' Sig0 = Sig0, y_all = y_all, c_all = OLCA_params$c_all)
#' 
#' # Then run MCMC sampling
#' MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, 
#' probit_params = probit_params, n_runs = 50, burn = 25, thin = 5,
#' K = K, p = p, d = d, n = n, q = q, w_all = w_all, x_mat = x_mat, 
#' y_all = y_all, V = V, alpha = alpha, eta = eta, Sig0 = Sig0, mu0 = mu0)
#' 
#' # Then run post-process relabeling
#' post_MCMC_out <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)
#'
#' # Then obtain posterior estimates
#' estimates <- get_estimates(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out,
#'                            n = n, p = p, V = V, y_all = y_all, x_mat = x_mat)
#' 
get_estimates <- function(MCMC_out, post_MCMC_out, n, p, V, y_all, x_mat) {
  
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
  pi_red <- post_MCMC_out$pi[, unique_classes] # Initialize pi for unique classes
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
  theta_red <- post_MCMC_out$theta[, , unique_classes, ]
  theta_red <- plyr::aaply(theta_red, c(1, 2, 3), function(x) x / sum(x)) # Re-normalize
  xi_red <- post_MCMC_out$xi[, unique_classes, , drop = FALSE]
  
  #============== Posterior median estimates ===================================
  pi_med <- apply(pi_red, 2, stats::median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), stats::median, na.rm = TRUE)
  theta_med <- plyr::aaply(theta_med, c(1, 2), function(x) x / sum(x))  # Re-normalize
  xi_med <- apply(xi_red, c(2, 3), stats::median, na.rm = TRUE)
  
  #============== Update c using unique classes and posterior estimates ========
  z_all <- MCMC_out$z_all_MCMC[M, ]
  c_all <- MCMC_out$c_all_MCMC[M, ]                      # Final class assignments
  pred_class_probs <- matrix(NA, nrow = n, ncol = K_red) # Posterior class membership probabilities
  log_cond_c <- matrix(NA, nrow = n, ncol = K_red)       # Individual log-likelihood for each class
  Phi_med_all_c <- stats::pnorm(V %*% t(xi_med))  # Outcome probabilities for all classes
  Phi_med <- numeric(n)                    # Initialize individual outcome probabilities
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    for (k in 1:K_red) {
      # Calculate theta component of individual log-likelihood assuming class k
      log_theta_comp_k <- 0
      for (j in 1:p) {
        log_theta_comp_k <- log_theta_comp_k + log(theta_med[j, k, x_mat[i, j]])
      }
      # Calculate and control extremes for probit component
      log_probit_part <- log(stats::dnorm(z_all[i], mean = V[i, ] %*% xi_med[k, ])) 
      if (log_probit_part == -Inf) {
        log_probit_part <- log(1e-16)
      }
      log_probit_comp_k <- log_probit_part + 
        log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
      # Individual log-likelihood for class k
      log_cond_c[i, k] <- log(pi_med[k]) + log_theta_comp_k + log_probit_comp_k
    }
    # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs[i, ] <- exp(log_cond_c[i, ] - matrixStats::logSumExp(log_cond_c[i, ]))
    # Update class assignment using the posterior probabilities
    c_all[i] <- LaplacesDemon::rcat(n = 1, p = pred_class_probs[i, ])
    # Calculate outcome probabilities P(Y=1|-) using updated class assignment
    Phi_med[i] <- Phi_med_all_c[i, c_all[i]]
  }
  
  #============== Update individual log-likelihood  ============================
  loglik_med <- numeric(n)  # Individual log-likelihood
  for (i in 1:n) {
    c_i <- c_all[i]
    # Calculate theta component of individual log-likelihood
    log_theta_comp <- 0
    for (j in 1:p) {
      log_theta_comp <- log_theta_comp + log(theta_med[j, c_i, x_mat[i, j]])
    }
    # Calculate individual log-likelihood using median estimates
    loglik_med[i] <- log(pi_med[c_i]) + log_theta_comp +
      log(stats::dnorm(z_all[i], mean = V[i, ] %*% xi_med[c_i, ])) + 
      log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
  }
  
  estimates <- list(K_red = K_red, pi_red = pi_red, theta_red = theta_red, 
                   xi_red = xi_red, pi_med = pi_med, theta_med = theta_med, 
                   xi_med = xi_med, Phi_med = Phi_med, c_all = c_all, 
                   pred_class_probs = pred_class_probs, loglik_med = loglik_med)
  return(estimates)
}