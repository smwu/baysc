#' Run MCMC to get posterior samples for WOLCA 
#'
#' @description
#' `run_MCMC_Rcpp_wolca` runs the Gibbs sampler MCMC algorithm using Rcpp to 
#' obtain posterior samples for the two-step unsupervised WOLCA model.
#'
#' @inheritParams run_MCMC_Rcpp
#' @inheritParams wolca
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
#' [wolca_svyglm()] [wolca()] 
#' @importFrom gtools permute
#' @export
#'
#' @examples
#'    
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
                                w_all, x_mat, alpha, eta, update = 10) {
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
    }
    
    # Print out progress 
    if (m %% update == 0) {
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
#' @param MCMC_out Output from `run_MCMC_Rcpp_wolca` containing:
#' \describe{
#'   \item{\code{pi_MCMC}}{Matrix of posterior samples for pi. (n_iter)xK}
#'   \item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xJxKxR}
#'   \item{\code{c_all_MCMC}}{Matrix of posterior samples for c_all. (n_iter)xn}
#' }
#' 
#' @details
#' First, `K_med`, the median number of classes with at least the `class_cutoff` 
#' proportion of individuals, is obtained over all MCMC iterations. Then, label 
#' switching is addressed through a relabeling procedure, where agglomerative 
#' clustering with Hamming distance is used to group individuals into `K_med` 
#' clusters and labels are re-assigned based on these clusters. Finally, 
#' parameter estimates are reordered using the relabeled classes so that 
#' posterior output can be averaged across MCMC iterations.
#' 
#' @return
#' Returns list `post_MCMC_out` containing:
#' \describe{
#'   \item{\code{K_med}}{Median, across iterations, of number of classes with at least `class_cutoff` percent of individuals}
#'   \item{\code{pi}}{Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)}
#'   \item{\code{theta}}{Array of reduced and relabeled posterior samples for theta. (n_iter)xJx(K_med)xR}
#'   \item{\code{dendrogram}}{Hierarchical clustering dendrogram used for relabeling}
#' }
#' 
#' @seealso [post_process()] [run_MCMC_Rcpp_wolca()] [get_estimates_wolca()] 
#' [wolca_svyglm()] [wolca()] 
#' @importFrom stats median as.dist hclust cutree
#' @importFrom e1071 hamming.distance
#' @export
#'
#' @examples
#'   
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
  # Modify classes if any classes are less than the cutoff percentage
  class_prop <- prop.table(table(red_c_all))
  if (any(class_prop < class_cutoff)) {
    # Get classes that are too small
    small <- which(class_prop < class_cutoff)
    # Group individuals into a larger number of classes 
    red_c_all_temp <- stats::cutree(dendrogram, k = K_med + length(small))
    red_c_all <- red_c_all_temp
    class_prop_temp <- prop.table(table(red_c_all_temp))
    # Get updated classes that are too small
    small_temp <- sort(which(class_prop_temp < class_cutoff))
    for (small_c in 1:length(small_temp)) {
      c_ind <- small_temp[small_c]
      class_small <- which(red_c_all_temp == c_ind)
      # Get nearest class
      inds <- 1:length(class_prop_temp)
      class_dist <- sapply(inds, function(x) 
        mean(distMat[class_small, which(red_c_all_temp == x)]))
      # Set small class distance to Inf
      class_dist[small_temp] <- Inf
      nearest <- which.min(class_dist[-c_ind])
      red_c_all[red_c_all_temp == c_ind] <- nearest
    }
    class_prop <- prop.table(table(red_c_all))
    K_med <- length(class_prop)
    # # Check class sizes
    # prop.table(table(red_c_all))
  }
  # Get unique reduced classes to aid relabeling
  unique_red_classes <- unique(red_c_all)
  
  # For each iteration, relabel new classes using the most common old class assignment
  relabel_red_classes <- matrix(NA, nrow = M, ncol = K_med)   # Apply new classes to each iteration
  for (k in 1:K_med) {
    red_class <- unique_red_classes[k]
    relabel_red_classes[, k] <- 
      apply(as.matrix(MCMC_out$c_all_MCMC[, red_c_all == red_class]), 1, get_mode)
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
#' @param MCMC_out Output from `run_MCMC_Rcpp_wolca` containing:
#' \describe{
#'   \item{\code{pi_MCMC}}{Matrix of posterior samples for pi. (n_iter)xK}
#'   \item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xJxKxR}
#'   \item{\code{c_all_MCMC}}{Matrix of posterior samples for c_all. (n_iter)xn}
#' }
#' @param post_MCMC_out output from `post_process_wolca` containing:
#' \describe{
#'   \item{\code{K_med}}{Median, across iterations, of number of classes with at least `class_cutoff` percent of individuals}
#'   \item{\code{pi}}{Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)}
#'   \item{\code{theta}}{Array of reduced and relabeled posterior samples for theta. (n_iter)xJx(K_med)xR}
#'   \item{\code{dendrogram}}{Hierarchical clustering dendrogram used for relabeling}
#' }
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
#' [wolca_svyglm()] [wolca()] 
#' @importFrom plyr aaply
#' @importFrom matrixStats logSumExp
#' @importFrom LaplacesDemon rcat
#' @importFrom stats median
#' @export
#'
#' @examples
#'    
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
  theta_red <- ifelse(theta_red < 1e-8, 1e-8, theta_red) # prevent underflow
  theta_red <- plyr::aaply(theta_red, c(1, 2, 3), function(x) x / sum(x), 
                           .drop = FALSE) # Re-normalize
  
  #============== Posterior median estimates ===================================
  pi_med <- apply(pi_red, 2, stats::median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), stats::median, na.rm = TRUE)
  theta_med <- ifelse(theta_med < 1e-8, 1e-8, theta_med) # prevent underflow
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



