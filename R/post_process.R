#' Post-process posterior samples to relabel and reduce classes
#'
#' @description
#' `post_process` conducts post-processing to cluster individuals into a 
#' reduced number of classes and reorder posterior parameter samples according 
#' to the reduced number of classes.
#'
#' @inheritParams run_MCMC_Rcpp
#' @inheritParams swolca
#' @param MCMC_out Output from `run_MCMC_Rcpp` containing `pi_MCMC`, 
#' `theta_MCMC`, `xi_MCMC`, `c_all_MCMC`, `z_all_MCMC`, and `loglik_MCMC`
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
#'   \item{\code{xi}}{Array of reduced and relabeled posterior samples for xi. (n_iter)x(K_med)xq}
#'   \item{\code{dendrogram}}{Hierarchical clustering dendrogram used for relabeling}
#' }
#' 
#' @seealso [run_MCMC_Rcpp()] [get_estimates()] [var_adjust()] [swolca()] [solca()]
#' @importFrom stats median as.dist hclust cutree
#' @importFrom e1071 hamming.distance
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
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
#' # Probit model only includes latent class
#' V <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
#' glm_form <- "~ 1"
#' # Obtain probit regression design matrix without class assignment
#' V <- model.matrix(as.formula(glm_form), data = V)
#' # Number of regression covariates excluding class assignment
#' q <- ncol(V)  
#' 
#' # Set hyperparameters
#' K <- 30
#' alpha <- rep(1, K) / K
#' eta <- matrix(0.01, nrow = J, ncol = R) 
#' for (j in 1:J) {
#'   eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
#' }
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
#' OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)
#' 
#' # Then initialize probit params 
#' probit_params <- init_probit(K = K, n = n, q = q, V = V, mu0 = mu0, 
#' Sig0 = Sig0, y_all = y_all, c_all = OLCA_params$c_all)
#' 
#' # Then run MCMC sampling
#' MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, 
#' probit_params = probit_params, n_runs = 50, burn = 25, thin = 5,
#' K = K, J = J, R = R, n = n, q = q, w_all = w_all, x_mat = x_mat, 
#' y_all = y_all, V = V, alpha = alpha, eta = eta, Sig0 = Sig0, mu0 = mu0)
#' 
#' # Then run post-process relabeling
#' post_MCMC_out <- post_process(MCMC_out = MCMC_out, J = J, R = R, q = q,
#' class_cutoff = 0.05)
#' # post_MCMC_out
#' # plot(post_MCMC_out$dendrogram)
#' 
post_process <- function(MCMC_out, J, R, q, class_cutoff) {
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
  xi <- array(NA, dim = c(M, K_med, q))
  for (m in 1:M) {
    iter_order <- relabel_red_classes[m, ]
    pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
    pi[m, ] <- pi_temp / sum(pi_temp)
    theta[m, , , ] <- MCMC_out$theta_MCMC[m, , iter_order, ]
    xi[m, , ] <- MCMC_out$xi_MCMC[m, iter_order, ]
  }
  
  post_MCMC_out <- list(K_med = K_med, pi = pi, theta = theta, xi = xi,
                        dendrogram = dendrogram)
  
  return(post_MCMC_out)  
  # plot(dendrogram, labels = FALSE)
  # rect.hclust(dendrogram, k = K_med)
}