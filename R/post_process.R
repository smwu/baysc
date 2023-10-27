#' Post-process posterior samples to relabel and reduce classes
#'
#' @description
#' `post_process` conducts post-processing to cluster individuals into a 
#' reduced number of classes and reorder posterior parameter samples according 
#' to the reduced number of classes.
#'
#' @inheritParams run_MCMC_Rcpp
#' @param MCMC_out Output from `run_MCMC_Rcpp` containing `pi_MCMC`, 
#' `theta_MCMC`, `xi_MCMC`, `c_all_MCMC`, `z_all_MCMC`, and `loglik_MCMC`
#' 
#' @return
#' Returns list `post_MCMC_out` containing:
#' \describe{
#'   \item{\code{K_med}}{Median, across iterations, of number of classes with at least 5 percent of individuals}
#'   \item{\code{pi}}{Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)}
#'   \item{\code{theta}}{Array of reduced and relabeled posterior samples for theta. (n_iter)xpx(K_med)xd}
#'   \item{\code{xi}}{Array of reduced and relabeled posterior samples for xi. (n_iter)x(K_med)xq}
#'   \item{\code{dendrogram}}{Hierarchical clustering dendrogram used for relabeling}
#' }
#' 
#' @importFrom stats median as.dist hclust cutree
#' @importFrom e1071 hamming.distance
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
#' # post_MCMC_out
#' # plot(post_MCMC_out$dendrogram)
#' 
post_process <- function(MCMC_out, p, d, q) {
  # Get median number of classes with >= 5% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(stats::median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
  
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
  theta <- array(NA, dim = c(M, p, K_med, d))
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