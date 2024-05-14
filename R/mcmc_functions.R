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
#' Finally, item level probabilities \eqn{\theta} are initialized by 
#' drawing from a Dirichlet distribution with hyperparameter \eqn{\eta}, 
#' independently for each exposure item and latent class.
#' 
#' @return
#' Returns list `OLCA_params` containing:
#' \describe{
#'   \item{\code{pi}}{Vector parameter pi for class membership probabilities. Kx1}
#'   \item{\code{theta}}{Array parameter theta for item level probabilities. JxKxR}
#'   \item{\code{c_all}}{Vector of random initial class assignments. nx1}
#' }
#' 
#' @seealso [init_probit()] [swolca()] [wolca()]
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



#' Initialize probit model
#'
#' @description
#' `init_probit` initializes priors and variables for the probit regression model 
#' given hyperparameters.
#'
#' @inheritParams init_OLCA
#' @param Q Number of regression covariates excluding class assignment
#' @param V Regression design matrix without class assignment. nxQ
#' @param mu0 List of K Qx1 vectors of mean hyperparameters for regression 
#' coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}. 
#' @param Sig0 List of K QxQ matrices of variance hyperparameters for regression 
#' coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}. 
#' @param y_all Vector of outcomes. nx1
#' @param c_all Vector of random initial class assignments. nx1
#' 
#' @details
#' First, regression coefficients \eqn{\xi} are initialized by drawing independently
#' for each latent class from a Multivariate Normal distribution with mean vector 
#' hyperparameter \eqn{\mu_0} drawn from a Normal(0,1) hyperprior and variance 
#' matrix hyperparameter \eqn{\Sigma_0} set to be a diagonal matrix with 
#' diagonal components drawn from InvGamma(shape=5/2, scale=2/5) distributions.
#' Then, the latent probit variable `z_all` is initialized by drawing from a 
#' Truncated Normal distribution with mean \eqn{V\xi} and variance 1 and 
#' truncation boundary 0, where negative values are drawn if the outcome is 0 and
#' positive values are drawn if the outcome is 1. 
#'
#' @return
#' Returns list `probit_params` containing:
#' \describe{
#'   \item{\code{xi}}{Matrix parameter xi for probit regression coefficients. KxQ}
#'   \item{\code{z_all}}{Vector of latent variables in the probit model. nx1}
#' }
#'
#' @seealso [init_OLCA()] [swolca()] 
#' 
#' @importFrom LaplacesDemon rmvn rinvgamma
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats qnorm rnorm
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' n <- dim(x_mat)[1]        # Number of individuals
#' J <- dim(x_mat)[2]        # Number of exposure items
#' R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
#'              function(x) length(unique(x)))  
#' R <- max(R_j)             # Maximum number of exposure categories across items
#' 
#' # Probit model only includes latent class
#' V_data <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
#' glm_form <- "~ 1"
#' # Obtain probit regression design matrix without class assignment
#' V <- model.matrix(as.formula(glm_form), data = V_data)
#' # Number of regression covariates excluding class assignment
#' Q <- ncol(V)  
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
#'   mu0[[k]] <- stats::rnorm(n = Q)
#'   # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
#'   # components and mean variance 2.5 for a weakly informative prior on xi
#'   Sig0[[k]] <- diag(LaplacesDemon::rinvgamma(n = Q, shape = 3.5, scale = 6.25), 
#'   nrow = Q, ncol = Q)
#' }
#' 
#' # First initialize OLCA params
#' OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)
#' 
#' # Initialize probit model
#' probit_params <- init_probit(K = K, n = n, Q = Q, V = V, mu0 = mu0, 
#' Sig0 = Sig0, y_all = y_all, c_all = OLCA_params$c_all)
#' # probit_params
#' 
init_probit <- function(K, n, Q, V, mu0, Sig0, y_all, c_all) {
  # Initialize variables
  xi <- matrix(NA, nrow = K, ncol = Q)
  z_all <- lin_pred <- numeric(n)
  
  # Prior for xi. Same prior for each class
  for (k in 1:K) {
    xi[k, ] <- LaplacesDemon::rmvn(n = 1, mu = mu0[[k]], Sigma = Sig0[[k]])
  }
  
  # Initialize probit model latent variable, z, following Albert and Chib (1993)
  # Linear predictor using covariate values and class assignment for each individual
  for (i in 1:n) {
    lin_pred[i] <- V[i, ] %*% xi[c_all[i], ]  # V[i]*xi[c_i]
  }
  # For cases, z_i ~ TruncNormal(low=0, high=Inf, mean=V*xi[c_i], var=1)
  z_all[y_all == 1] <- truncnorm::rtruncnorm(n = 1, a = 0, b = Inf, 
                                             mean = lin_pred[y_all == 1], sd = 1)
  # For controls, z_i ~ TruncNormal(low=-Inf, high=0, mean=V*xi[c_i], var=1)
  z_all[y_all == 0] <- truncnorm::rtruncnorm(n = 1, a = -Inf, b = 0, 
                                             mean = lin_pred[y_all == 0], sd = 1)
  # Control extreme values
  z_all[z_all == Inf] <- stats::qnorm(1 - (1e-10))
  z_all[z_all == -Inf] <- stats::qnorm(1e-10)
  
  probit_params <- list(xi = xi, z_all = z_all)
  return(probit_params)
}



#' Run MCMC to get posterior samples
#'
#' @description
#' `run_MCMC_Rcpp` runs the Gibbs sampler MCMC algorithm using Rcpp to obtain 
#' posterior samples.
#'
#' @inheritParams init_OLCA
#' @inheritParams init_probit
#' @inheritParams swolca
#' @param OLCA_params Output list from [init_OLCA()] containing:
#' \describe{
#'   \item{\code{pi}}{Vector parameter pi for class membership probabilities. Kx1}
#'   \item{\code{theta}}{Array parameter theta for item level probabilities. JxKxR}
#'   \item{\code{c_all}}{Vector of random initial class assignments. nx1}
#' }
#' @param probit_params Output list from [init_probit()] containing:
#' \describe{
#'   \item{\code{xi}}{Matrix parameter xi for probit regression coefficients. KxQ}
#'   \item{\code{z_all}}{Vector of latent variables in the probit model. nx1}
#' }
#' @param w_all Weights normalized to sum to n. nx1
#' 
#' @details
#' A Gibbs sampler updates the parameters and variables in the following order:
#' \eqn{\pi}, `c_all`, \eqn{\theta}, \eqn{\xi}, `z_all`. Class assignments
#' are permuted every 10 iterations to encourage mixing, according to a random
#' permutation sampler (Fruhwirth-Schnatter, 2001).
#' 
#' @return
#' Returns list `MCMC_out` containing:
#' \describe{
#'   \item{\code{pi_MCMC}}{Matrix of posterior samples for pi. (n_iter)xK}
#'   \item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xJxKxR}
#'   \item{\code{xi_MCMC}}{Array of posterior samples for xi. (n_iter)xKxQ}
#'   \item{\code{c_all_MCMC}}{Matrix of posterior samples for c_all. (n_iter)xn}
#'   \item{\code{z_all_MCMC}}{Matrix of posterior samples for z_all. (n_iter)xn}
#'   \item{\code{loglik_MCMC}}{Vector of posterior samples for log-likelihood. (n_iter)x1}
#' }
#'
#' @seealso [post_process()] [get_estimates()] [swolca_var_adjust()] [swolca()] 
#' [run_MCMC_Rcpp_wolca()]
#' @importFrom gtools permute
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats rnorm
#' @export
#' 
#' @references 
#' Fruhwirth-Schnatter, S. (2001). Markov chain monte carlo estimation of 
#' classical and dynamic switching and mixture models. Journal of the American 
#' Statistical Association 96, 194–209.
#' 
#' @examples
#'    
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
#' V_data <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
#' glm_form <- "~ 1"
#' # Obtain probit regression design matrix without class assignment
#' V <- model.matrix(as.formula(glm_form), data = V_data)
#' # Number of regression covariates excluding class assignment
#' Q <- ncol(V)  
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
#'   mu0[[k]] <- stats::rnorm(n = Q)
#'   # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
#'   # components and mean variance 2.5 for a weakly informative prior on xi
#'   Sig0[[k]] <- diag(LaplacesDemon::rinvgamma(n = Q, shape = 3.5, scale = 6.25), 
#'   nrow = Q, ncol = Q)
#' }
#' 
#' # First initialize OLCA params
#' OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)
#' 
#' # Then initialize probit params 
#' probit_params <- init_probit(K = K, n = n, Q = Q, V = V, mu0 = mu0, 
#' Sig0 = Sig0, y_all = y_all, c_all = OLCA_params$c_all)
#' 
#' # Then run MCMC sampling
#' MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, 
#' probit_params = probit_params, n_runs = 50, burn = 25, thin = 5,
#' K = K, J = J, R = R, n = n, Q = Q, w_all = w_all, x_mat = x_mat, 
#' y_all = y_all, V = V, alpha = alpha, eta = eta, Sig0 = Sig0, mu0 = mu0)
#' # MCMC_out
#' 
run_MCMC_Rcpp <- function(OLCA_params, probit_params, n_runs, burn, thin, K, J, 
                          R, n, Q, w_all, x_mat, y_all, V, alpha, eta, mu0, Sig0, 
                          update = 10) {
  # Number of MCMC iterations to store
  n_storage <- floor(n_runs / thin) 
  
  # Initialize variables
  pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
  theta_MCMC <- array(NA, dim = c(n_storage, J, K, R))
  xi_MCMC <- array(NA, dim = c(n_storage, K, Q))
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
    update_c(c_all = c_all, n = n, K = K, J = J, theta = theta, 
             x_mat = x_mat, pi = pi, z_all = z_all, V = V, xi = xi, 
             y_all = y_all)
    update_theta(theta = theta, J = J, K = K, R = R, eta = eta, 
                 w_all = w_all, c_all = c_all, x_mat = x_mat)
    xi <- update_xi(xi = xi, n = n, K = K, w_all = w_all, c_all = c_all,
                    z_all = z_all, V = V, mu0 = mu0, Sig0 = Sig0)
    z_all <- update_z(z_all = z_all, n = n, V = V, xi = xi, c_all = c_all,
                      y_all = y_all)
    update_loglik(loglik = loglik, n = n, J = J, c_all = c_all, 
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
      new_order <- gtools::permute(1:K)      # New permuted labels
      new_c_all <- numeric(K)        # Class assignments with new labels
      for (k in 1:K) {
        new_c_all[c_all == k] <- new_order[k]
      }
      c_all <- new_c_all                # Relabel class assignments
      pi <- pi[new_order]               # Relabel class probabilities
      theta <- theta[, new_order, , drop = FALSE]     # Relabel item category probabilities
      xi <- xi[new_order, , drop = FALSE]  # Relabel probit coefficients
    }
    
    # Print out progress 
    if (m %% update == 0) {
      print(paste0("Iteration ", m, " completed!"))
    }
  }
  
  # Discard burn-in
  warmup <- floor(burn / thin)
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



#' Post-process posterior samples to relabel and reduce classes
#'
#' @description
#' `post_process` conducts post-processing to cluster individuals into a 
#' reduced number of classes and reorder posterior parameter samples according 
#' to the reduced number of classes.
#'
#' @inheritParams run_MCMC_Rcpp
#' @inheritParams swolca
#' @param MCMC_out Output from [run_MCMC_Rcpp()] containing:
#' \describe{
#'   \item{\code{pi_MCMC}}{Matrix of posterior samples for pi. (n_iter)xK}
#'   \item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xJxKxR}
#'   \item{\code{xi_MCMC}}{Array of posterior samples for xi. (n_iter)xKxQ}
#'   \item{\code{c_all_MCMC}}{Matrix of posterior samples for c_all. (n_iter)xn}
#'   \item{\code{z_all_MCMC}}{Matrix of posterior samples for z_all. (n_iter)xn}
#'   \item{\code{loglik_MCMC}}{Vector of posterior samples for log-likelihood. (n_iter)x1}
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
#'   \item{\code{xi}}{Array of reduced and relabeled posterior samples for xi. (n_iter)x(K_med)xQ}
#'   \item{\code{dendrogram}}{Hierarchical clustering dendrogram used for relabeling}
#' }
#' 
#' @seealso [run_MCMC_Rcpp()] [get_estimates()] [swolca_var_adjust()] [swolca()] 
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
#' V_data <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
#' glm_form <- "~ 1"
#' # Obtain probit regression design matrix without class assignment
#' V <- model.matrix(as.formula(glm_form), data = V_data)
#' # Number of regression covariates excluding class assignment
#' Q <- ncol(V)  
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
#'   mu0[[k]] <- stats::rnorm(n = Q)
#'   # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
#'   # components and mean variance 2.5 for a weakly informative prior on xi
#'   Sig0[[k]] <- diag(LaplacesDemon::rinvgamma(n = Q, shape = 3.5, scale = 6.25), 
#'   nrow = Q, ncol = Q)
#' }
#' 
#' # First initialize OLCA params
#' OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)
#' 
#' # Then initialize probit params 
#' probit_params <- init_probit(K = K, n = n, Q = Q, V = V, mu0 = mu0, 
#' Sig0 = Sig0, y_all = y_all, c_all = OLCA_params$c_all)
#' 
#' # Then run MCMC sampling
#' MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, 
#' probit_params = probit_params, n_runs = 50, burn = 25, thin = 5,
#' K = K, J = J, R = R, n = n, Q = Q, w_all = w_all, x_mat = x_mat, 
#' y_all = y_all, V = V, alpha = alpha, eta = eta, Sig0 = Sig0, mu0 = mu0)
#' 
#' # Then run post-process relabeling
#' post_MCMC_out <- post_process(MCMC_out = MCMC_out, J = J, R = R, Q = Q,
#' class_cutoff = 0.05)
#' # post_MCMC_out
#' # plot(post_MCMC_out$dendrogram)
#' 
post_process <- function(MCMC_out, J, R, Q, class_cutoff) {
  # Get median number of classes with >= cutoff% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(stats::median(rowSums(MCMC_out$pi_MCMC >= class_cutoff)))
  
  # Cluster individuals into reduced number of classes using agglomerative clustering
  # Calculate pairwise distance matrix using Hamming distance: number of 
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
  xi <- array(NA, dim = c(M, K_med, Q))
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



#' Get posterior estimates
#'
#' @description
#' `get_estimates` obtains posterior parameter samples and estimates prior to
#' variance adjustment
#'
#' @inheritParams run_MCMC_Rcpp
#' @param MCMC_out Output from [run_MCMC_Rcpp()] containing:
#' \describe{
#'   \item{\code{pi_MCMC}}{Matrix of posterior samples for pi. (n_iter)xK}
#'   \item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xJxKxR}
#'   \item{\code{xi_MCMC}}{Array of posterior samples for xi. (n_iter)xKxQ}
#'   \item{\code{c_all_MCMC}}{Matrix of posterior samples for c_all. (n_iter)xn}
#'   \item{\code{z_all_MCMC}}{Matrix of posterior samples for z_all. (n_iter)xn}
#'   \item{\code{loglik_MCMC}}{Vector of posterior samples for log-likelihood. (n_iter)x1}
#' }
#' @param post_MCMC_out output from [post_process()] containing:
#' \describe{
#'   \item{\code{K_med}}{Median, across iterations, of number of classes with at least `class_cutoff` percent of individuals}
#'   \item{\code{pi}}{Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)}
#'   \item{\code{theta}}{Array of reduced and relabeled posterior samples for theta. (n_iter)xJx(K_med)xR}
#'   \item{\code{xi}}{Array of reduced and relabeled posterior samples for xi. (n_iter)x(K_med)xQ}
#'   \item{\code{dendrogram}}{Hierarchical clustering dendrogram used for relabeling}
#' }
#' 
#' @details
#' First, duplicate classes that have the same modal exposure categories
#' for all items are combined to obtain the number of unique classes, `K_red`. 
#' Parameters are then renormalized for the unique classes and posterior median 
#' estimates are computed across MCMC iterations. Using these median estimates,
#' class assignments `c_all`, the regression mean, and the individual 
#' log-likelihood are derived. 
#' 
#' @return
#' Returns list `estimates` containing:
#' \describe{
#'   \item{\code{K_red}}{Number of unique classes}
#'   \item{\code{pi_red}}{Matrix of final posterior samples for pi. Mx(K_red), 
#'   where M is the number of MCMC iterations after burn-in and thinning.}
#'   \item{\code{theta_red}}{Array of final posterior samples for theta. MxJx(K_red)xR}
#'   \item{\code{xi_red}}{Array of final posterior samples for xi. Mx(K_red)xQ}
#'   \item{\code{pi_med}}{Vector of posterior median estimates for pi. (K_red)x1}
#'   \item{\code{theta_med}}{Array of posterior median estimates for theta. Jx(K_red)xR}
#'   \item{\code{xi_med}}{Matrix of posterior median estimates for xi. (K_red)xQ}
#'   \item{\code{Phi_med}}{Vector of final individual outcome probabilities. nx1}
#'   \item{\code{c_all}}{Vector of final individual class assignments. nx1}
#'   \item{\code{pred_class_probs}}{Matrix of individual posterior class probabilities. nx(K_red)}
#'   \item{\code{loglik_med}}{Vector of final indiviudal log-likehoods. nx1} 
#' }
#'
#' @seealso [run_MCMC_Rcpp()] [post_process()] [swolca_var_adjust()] [swolca()] 
#' @importFrom plyr aaply
#' @importFrom matrixStats logSumExp
#' @importFrom LaplacesDemon rcat
#' @importFrom stats dnorm median pnorm
#' @export
#' 
#' @references 
#' Williams, M. R. and Savitsky, T. D. (2021). Uncertainty estimation for 
#' pseudo-bayesian inference under complex sampling. International Statistical 
#' Review 89, 72–107.
#'
#' @examples
#'    
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
#' V_data <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
#' glm_form <- "~ 1"
#' # Obtain probit regression design matrix without class assignment
#' V <- model.matrix(as.formula(glm_form), data = V_data)
#' # Number of regression covariates excluding class assignment
#' Q <- ncol(V)  
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
#'   mu0[[k]] <- stats::rnorm(n = Q)
#'   # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
#'   # components and mean variance 2.5 for a weakly informative prior on xi
#'   Sig0[[k]] <- diag(LaplacesDemon::rinvgamma(n = Q, shape = 3.5, scale = 6.25), 
#'   nrow = Q, ncol = Q)
#' }
#' 
#' # First initialize OLCA params
#' OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)
#' 
#' # Then initialize probit params 
#' probit_params <- init_probit(K = K, n = n, Q = Q, V = V, mu0 = mu0, 
#' Sig0 = Sig0, y_all = y_all, c_all = OLCA_params$c_all)
#' 
#' # Then run MCMC sampling
#' MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, 
#' probit_params = probit_params, n_runs = 50, burn = 25, thin = 5,
#' K = K, J = J, R = R, n = n, Q = Q, w_all = w_all, x_mat = x_mat, 
#' y_all = y_all, V = V, alpha = alpha, eta = eta, Sig0 = Sig0, mu0 = mu0)
#' 
#' # Then run post-process relabeling
#' post_MCMC_out <- post_process(MCMC_out = MCMC_out, J = J, R = R, Q = Q,
#' class_cutoff = 0.05)
#'
#' # Then obtain posterior estimates
#' estimates <- get_estimates(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out,
#'                            n = n, J = J, V = V, y_all = y_all, x_mat = x_mat)
#' 
get_estimates <- function(MCMC_out, post_MCMC_out, n, J, V, y_all, x_mat) {
  
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
  theta_red <- post_MCMC_out$theta[, , unique_classes, ]
  theta_red <- ifelse(theta_red < 1e-8, 1e-8, theta_red) # prevent underflow
  theta_red <- plyr::aaply(theta_red, c(1, 2, 3), function(x) x / sum(x),
                           .drop = FALSE) # Re-normalize
  xi_red <- post_MCMC_out$xi[, unique_classes, , drop = FALSE]
  
  #============== Posterior median estimates ===================================
  pi_med <- apply(pi_red, 2, stats::median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), stats::median, na.rm = TRUE)
  theta_med <- ifelse(theta_med < 1e-8, 1e-8, theta_med) # prevent underflow
  theta_med <- plyr::aaply(theta_med, c(1, 2), function(x) x / sum(x),
                           .drop = FALSE)  # Re-normalize
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
      for (j in 1:J) {
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
    for (j in 1:J) {
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