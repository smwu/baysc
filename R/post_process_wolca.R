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