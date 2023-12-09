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
#'   \item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xJxKxR}
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
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
#' cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
#' sampling_wt <- data_vars$sample_wt
#' 
#' # Obtain dimensions
#' n <- dim(x_mat)[1]        # Number of individuals
#' J <- dim(x_mat)[2]        # Number of exposure items
#' R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
#'              function(x) length(unique(x)))  
#' R <- max(R_j)             # Maximum number of exposure categories across items
#' # Obtain normalized weights
#' kappa <- sum(sampling_wt) / n   # Weights norm. constant
#' w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
#' 
#' # Set hyperparameters
#' K <- 30
#' alpha <- rep(1, K) / K
#' eta <- matrix(0.01, nrow = J, ncol = R) 
#' for (j in 1:J) {
#'   eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
#' }
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



#' Post-process posterior samples to relabel and reduce classes for WOLCA 
#'
#' @description
#' `post_process_wolca` conducts post-processing to cluster individuals into a 
#' reduced number of classes and reorder posterior parameter samples according 
#' to the reduced number of classes for the unsupervised WOLCA model.
#'
#' @inheritParams run_MCMC_Rcpp_wolca
#' @inheritParams wolca
#' @param MCMC_out Output from `run_MCMC_Rcpp_wolca` containing `pi_MCMC`, 
#' `theta_MCMC`, and `c_all_MCMC`
#' 
#' @details
#' First, `K_med`, the median number of classes with at least the `class_cutoff` 
#' proportion of individuals is obtained over all MCMC iterations. Then, label 
#' switching is addressed through a relabeling procedure, where agglomerative 
#' clustering with Hamming distance is used to group individuals into `K_med` 
#' clusters and labels are re-assigned based on these clusters. Finally, 
#' parameter estimates are reordered using the relabeled classes so that 
#' posterior output can be averaged across MCMC iterations.
#' 
#' @return
#' Returns list `post_MCMC_out` containing:
#' \describe{
#'   \item{\code{K_med}}{Median, across iterations, of number of classes with at least 5 percent of individuals}
#'   \item{\code{pi}}{Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)}
#'   \item{\code{theta}}{Array of reduced and relabeled posterior samples for theta. (n_iter)xJx(K_med)xR}
#'   \item{\code{dendrogram}}{Hierarchical clustering dendrogram used for relabeling}
#' }
#' 
#' @seealso [post_process()] [run_MCMC_Rcpp_wolca()] [get_estimates_wolca()] 
#' [fit_probit_wolca()] [wolca()] 
#' @importFrom stats median as.dist hclust cutree
#' @importFrom e1071 hamming.distance
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
#' cluster_id <- data_vars$cluster_id  # Cluster indicators, nx1
#' sampling_wt <- data_vars$sample_wt
#' 
#' # Obtain dimensions
#' n <- dim(x_mat)[1]        # Number of individuals
#' J <- dim(x_mat)[2]        # Number of exposure items
#' R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
#'              function(x) length(unique(x)))  
#' R <- max(R_j)             # Maximum number of exposure categories across items 
#' # Obtain normalized weights
#' kappa <- sum(sampling_wt) / n   # Weights norm. constant
#' w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
#' 
#' # Set hyperparameters for fixed sampler
#' K <- 3
#' alpha <- rep(1, K) / K
#' eta <- matrix(0.01, nrow = J, ncol = R) 
#' for (j in 1:J) {
#'   eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
#' }
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
#' # post_MCMC_out
#' # plot(post_MCMC_out$dendrogram)
#' 
post_process_wolca <- function(MCMC_out, J, R, class_cutoff) {
  # Get median number of classes with >= cutoff% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(stats::median(rowSums(MCMC_out$pi_MCMC >= class_cutoff)))
  
  # Cluster individuals into reduced number of classes using agglomerative clustering
  # Calculate pairwise distance matrix using Hamming distance: proportion of
  # iterations where two individuals have differing class assignments
  distMat <- e1071::hamming.distance(t(MCMC_out$c_all_MCMC))
  # Hierarchical clustering dendrogram
  dendrogram <- stats::hclust(stats::as.dist(distMat), method = "complete") 
  # Group individuals into K_med classes
  red_c_all <- stats::cutree(dendrogram, k = K_med)                  
  # For each iteration, relabel new classes using the most common old class assignment
  relabel_red_classes <- matrix(NA, nrow = M, ncol = K_med)   # Apply new classes to each iteration
  for (k in 1:K_med) {
    relabel_red_classes[, k] <- apply(as.matrix(MCMC_out$c_all_MCMC[, red_c_all == k]),
                                      1, get_mode)
  }
  
  # Reduce and reorder parameter estimates using new classes
  pi <- matrix(NA, nrow = M, ncol = K_med)
  theta <- array(NA, dim = c(M, J, K_med, R))
  for (m in 1:M) {
    iter_order <- relabel_red_classes[m, ]
    pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
    pi[m, ] <- pi_temp / sum(pi_temp)
    theta[m, , , ] <- MCMC_out$theta_MCMC[m, , iter_order, ]
  }
  
  post_MCMC_out <- list(K_med = K_med, pi = pi, theta = theta,
                        dendrogram = dendrogram)
  
  return(post_MCMC_out)
  # plot(dendrogram, labels = FALSE)
  # rect.hclust(dendrogram, k = K_med)
}



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
#'   \item{\code{theta_red}}{Array of final posterior samples for theta. MxJx(K_red)xR}
#'   \item{\code{pi_med}}{Vector of posterior median estimates for pi. (K_red)x1}
#'   \item{\code{theta_med}}{Array of posterior median estimates for theta. Jx(K_red)xR}
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
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
#' cluster_id <- data_vars$cluster_id  # Cluster indicators, nx1
#' sampling_wt <- data_vars$sample_wt
#' 
#' # Obtain dimensions
#' n <- dim(x_mat)[1]        # Number of individuals
#' J <- dim(x_mat)[2]        # Number of exposure items
#' R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
#'              function(x) length(unique(x)))  
#' R <- max(R_j)             # Maximum number of exposure categories across items
#' # Obtain normalized weights
#' kappa <- sum(sampling_wt) / n   # Weights norm. constant
#' w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
#' 
#' # Set hyperparameters for fixed sampler
#' K <- 3
#' alpha <- rep(1, K) / K
#' eta <- matrix(0.01, nrow = J, ncol = R) 
#' for (j in 1:J) {
#'   eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
#' }
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



#' Fit probit model for WOLCA
#'
#' @description
#' `fit_probit_wolca` uses `svyglm` to fit a survey-weight probit model in a 
#' two-step model where the first step derived latent classes using an 
#' unsupervised WOLCA 
#'
#' @inheritParams wolca
#' @param estimates Output from `get_estimates_wolca()` containing `K_red`, 
#' `pi_red`, `theta_red`, `pi_med`, `theta_med`, `c_all`, `pred_class_probs`
#' @param w_all Weights normalized to sum to n. nx1
#' @param q Number of regression covariates excluding class assignment
#' 
#' @details
#' Specifies survey design and fits a survey-weighted probit regression model
#' according to the formula specified in `glm_form`. Regression coefficients and 
#' their confidence intervals are obtained from the `svyglm()` output. If the 
#' residual degrees of freedom is less than 1, a Wald confidence interval is 
#' manually calculated using a t-distribution with degrees of freedom from the 
#' survey design. The point and interval estimates are then converted into the 
#' factor reference coding format to match the output from `swolca()` and `solca()`. 
#' 
#' @return
#' Returns updated list `estimates` containing the following additional objects:
#' \describe{
#'   \item{\code{xi_est}}{Matrix of estimates for xi. (K_red)xq}
#'   \item{\code{xi_est_lb}}{Matrix of confidence interval lower bound estimates for xi. (K_red)xq}
#'   \item{\code{xi_est_ub}}{Matrix of confidence interval upper bound estimates for xi. (K_red)xq}
#'   \item{\code{fit}}{`svyglm` class object with output from the `svyglm` regression model}
#' }
#'
#' @seealso [run_MCMC_Rcpp_wolca()] [post_process_wolca()] 
#' [get_estimates_wolca()] [wolca()] 
#' @importFrom survey svydesign svyglm degf
#' @importFrom stats confint as.formula quasibinomial terms as.formula
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' cluster_id <- data_vars$cluster_id  # Cluster indicators, nx1
#' stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
#' sampling_wt <- data_vars$sample_wt
#' 
#' # Obtain dimensions
#' n <- dim(x_mat)[1]        # Number of individuals
#' J <- dim(x_mat)[2]        # Number of exposure items
#' R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
#'              function(x) length(unique(x)))  
#' R <- max(R_j)             # Maximum number of exposure categories across items
#' # Obtain normalized weights
#' kappa <- sum(sampling_wt) / n   # Weights norm. constant
#' w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
#' 
#' # Set hyperparameters for fixed sampler
#' K <- 3
#' alpha <- rep(1, K) / K
#' eta <- matrix(0.01, nrow = J, ncol = R) 
#' for (j in 1:J) {
#'   eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
#' }
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
#' # Then obtain posterior estimates for WOLCA
#' estimates <- get_estimates_wolca(MCMC_out = MCMC_out, 
#' post_MCMC_out = post_MCMC_out, n = n, J = J, x_mat = x_mat)
#' 
#' # Probit model only includes latent class
#' V_data <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
#' # Survey-weighted regression formula
#' glm_form <- "~ 1"
#' # Obtain probit regression design matrix without class assignment
#' V <- model.matrix(as.formula(glm_form), data = V_data)
#' # Number of regression covariates excluding class assignment
#' q <- ncol(V)  
#' 
#' # Finally run weighted probit regression model
#' estimates <- fit_probit_wolca(estimates = estimates, glm_form = glm_form, 
#' stratum_id = stratum_id, cluster_id = cluster_id, x_mat = x_mat, 
#' y_all = y_all, w_all = w_all, V_data = V_data, q = q)
#' 
fit_probit_wolca <- function(estimates, glm_form, stratum_id, cluster_id, 
                             x_mat, y_all, w_all, ci_level = 0.95, V_data, q) {
  
  # Create survey design
  if (!is.null(stratum_id)) {  # Include stratifying variable
    # Survey data frame for specifying survey design
    svy_data <- data.frame(stratum_id = factor(stratum_id), 
                           cluster_id = factor(cluster_id),
                           x_mat = x_mat, y_all = y_all, w_all = w_all)
    # Add latent class assignment variable to survey data
    svy_data$c_all <- factor(estimates$c_all)
    # Add additional covariates
    svy_data <- cbind(svy_data, V_data)
    # Specify survey design
    svydes <- survey::svydesign(ids = ~cluster_id, strata = ~stratum_id, 
                                weights = ~w_all, data = svy_data)
  } else { # No stratifying variable
    # Survey data frame for specifying survey design
    svy_data <- data.frame(cluster_id = cluster_id, x_mat = x_mat, 
                           y_all = y_all, w_all = w_all)
    # Add latent class assignment variable to survey data
    svy_data$c_all <- factor(estimates$c_all)
    # Add additional covariates
    svy_data <- cbind(svy_data, V_data)
    # Specify survey design
    svydes <- survey::svydesign(ids = ~cluster_id, weights = ~w_all, 
                                data = svy_data)    
  }
  
  # Add outcome and latent class main and interaction terms to formula
  terms <- labels(stats::terms(stats::as.formula(glm_form)))
  if (length(terms) > 1) {
    full_glm_form <- paste0("y_all ~ ", 
                            paste0("c_all * ", terms, collapse = " + ")) 
  } else {
    full_glm_form <- paste0("y_all ~ c_all") 
  }
  
  # If only one latent class, cannot have latent class as a covariate
  if (length(levels(svy_data$c_all)) == 1) {
    stop("Only one latent class found. Cannot use latent class as a covariate")
  }
  
  # Fit probit model according to specified formula
  fit <- survey::svyglm(formula = stats::as.formula(full_glm_form), 
                        design = svydes, 
                        family = stats::quasibinomial(link = "probit"))
  # Obtain coefficients and confidence interval
  coefs <- fit$coefficients
  ci <- stats::confint(fit, level = ci_level)
  # If zero/negative residual df, manually calculate the Wald confidence interval 
  # using a t-distribution with degrees of freedom from the survey design. 
  # Best if no cluster-level covariates in the regression model
  if (all(is.na(ci))) {
    ci <- manual_CI(model_object = fit, svy_df = survey::degf(svydes), 
                    ci = ci_level)[, -1]
  }
  
  # Convert format to match SWOLCA and SOLCA
  xi_list <- convert_ref_to_mix(K = estimates$K_red, q = q, est_beta = coefs,
                                ci_beta = ci)
  
  # Return output with probit regression estimates
  estimates$xi_est <- xi_list$est_xi
  estimates$xi_est_lb <- xi_list$est_xi_lb
  estimates$xi_est_ub <- xi_list$est_xi_ub
  estimates$fit <- fit
  
  return(estimates)
}
