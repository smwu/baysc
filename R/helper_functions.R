#===================================================
## Helper functions for WSOLCA, SOLCA, and WOLCA
## Programmer: SM Wu   
## Data: Simulations and application   
#===================================================

# 'init_OLCA' initializes priors and variables for the OLCA model given hyperparameters
# Inputs:
#   alpha: Vector of hyperparameters for pi. Kx1
#   eta: Vector of hyperparameters for theta. dx1
#   n: Number of individuals
#   K: Number of classes
#   p: Number of exposure items
#   d: Number of exposure categories
# Outputs: `OLCA_params` list with the following items:
#   pi: Vector parameter pi for class membership probabilities. Kx1
#   theta: Array parameter theta for item category probabilities. pxKxd
#   c_all: Vector of random initial class assignments. nx1
init_OLCA <- function(alpha, eta, n, K, p, d) {
  # Prior for pi
  pi <- c(rdirichlet(n = 1, alpha = alpha))
  
  # Initialize class assignment, c, for individuals
  c_all <- rcat(n = n, p = pi)
  
  # Prior for theta
  theta <- array(0, dim = c(p, K, d))
  for (j in 1:p) {
    for (k in 1:K) {
      theta[j, k, ] <- c(rdirichlet(n = 1, alpha = eta))
    }
  }
  
  OLCA_params <- list(pi = pi, theta = theta, c_all = c_all)
  return(OLCA_params)
}

# `init_probit` initializes priors and variables for the probit regression model
# given hyperparameters 
# Inputs:
#   mu0: List of vectors of mean hyperparameters for xi. List of K qx1 vectors
#   Sig0: List of matrices of variance hyperparameters for xi. List of K qxq matrices
#   K: Number of classes
#   q: Number of regression covariates excluding class assignment
#   n: Number of individuals
#   V: Regression design matrix without class assignment. nxq
#   y_all: Vector of outcomes. nx1
#   c_all: Vector of random initial class assignments. nx1
# Outputs: `probit_params` list with the following items:
#   xi: Matrix parameter xi for probit regression coefficients. Kxq
#   z_all: Vector of latent variables in the probit model. nx1
init_probit <- function(mu0, Sig0, K, q, n, V, y_all, c_all) {
  # Initialize variables
  xi <- matrix(NA, nrow = K, ncol = q)
  z_all <- lin_pred <- numeric(n)
  
  # Prior for xi. Same prior for each class
  for (k in 1:K) {
    xi[k, ] <- rmvn(n = 1, mu = mu0[[k]], Sigma = Sig0[[k]])
  }
  
  # Initialize probit model latent variable, z, following Albert and Chib (1993)
  # Linear predictor using covariate values and class assignment for each individual
  for (i in 1:n) {
    lin_pred[i] <- V[i, ] %*% xi[c_all[i], ]  # V[i]*xi[c_i]
  }
  # For cases, z_i ~ TruncNormal(low=0, high=Inf, mean=V*xi[c_i], var=1)
  z_all[y_all == 1] <- rtruncnorm(n = 1, a = 0, b = Inf, 
                                  mean = lin_pred[y_all == 1], sd = 1)
  # For controls, z_i ~ TruncNormal(low=-Inf, high=0, mean=V*xi[c_i], var=1)
  z_all[y_all == 0] <- rtruncnorm(n = 1, a = -Inf, b = 0, 
                                  mean = lin_pred[y_all == 0], sd = 1)
  # Control extreme values
  z_all[z_all == Inf] <- qnorm(1 - (1e-10))
  z_all[z_all == -Inf] <- qnorm(1e-10)
  
  probit_params <- list(xi = xi, z_all = z_all)
  return(probit_params)
}

update_xi_R <- function(n, w_all, V, K, c_all, Sig0, mu0, z_all, xi) {
  W_tilde <- sparseMatrix(i = 1:n, j = 1:n, x = w_all)  # Sparse diagonal normalized weight matrix
  V_sparse <- as(V, "sparseMatrix")                     # Sparse design matrix without class assignments
  for (k in 1:K) {
    # Sparse diagonal matrix subsetting to obs in class k
    C_k <- sparseMatrix(i = 1:n, j = 1:n, x = (c_all == k))
    
    # Draw xi from conditional posterior distribution
    Sig_post <- solve(Sig0[[k]]) + as.matrix(t(V) %*% C_k %*% W_tilde %*% V)
    mu_right <- (solve(Sig0[[k]]) %*% mu0[[k]]) +
      (as.matrix(t(V) %*% C_k %*% W_tilde %*% z_all))
    mu_post <- c(solve(Sig_post) %*% mu_right)
    xi[k, ] <- rmvn(n = 1, mu = mu_post, Sigma = solve(Sig_post))
  }
  return(xi)
}

update_z_R <- function(lin_pred, n, V, xi, c_all, z_all, y_all) {
  # Linear predictor using covariate values and class assignment for each individual
  for (i in 1:n) {
    lin_pred[i] <- V[i, ] %*% xi[c_all[i], ]  # V[i]*xi[c_i]
  }
  # Probit model latent variable z, following Albert and Chib (1993)
  # For cases, z_i ~ TruncNormal(low=0, high=Inf, mean=V*xi[c_i], var=1)
  z_all[y_all == 1] <- rtruncnorm(n = 1, a = 0, b = Inf,
                                  mean = lin_pred[y_all == 1], sd = 1)
  # For controls, z_i ~ TruncNormal(low=-Inf, high=0, mean=V*xi[c_i], var=1)
  z_all[y_all == 0] <- rtruncnorm(n = 1, a = -Inf, b = 0,
                                  mean = lin_pred[y_all == 0], sd = 1)
  # Control extreme values
  z_all[z_all == Inf] <- qnorm(1 - (1e-10))
  z_all[z_all == -Inf] <- qnorm(1e-10)
  return(z_all)
}



# `run_MCMC_Rcpp` runs the Gibbs sampler MCMC algorithm to obtain posterior samples
# Inputs:
#   OLCA_params: output from 'init_OLCA' containing 'pi', 'c_all', 'theta'
#   probit_params: output from 'init_probit' containing 'xi', 'z_all'
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
#   K: Number of classes
#   p: Number of exposure items
#   d: Number of exposure categories
#   n: Number of individuals
#   q: Number of regression covariates excluding class assignment
#   w_all: Weights normalized to sum to n. nx1
#   x_mat: Categorical exposure matrix. nxp
#   y_all: Vector of outcomes. nx1
#   V: Regression design matrix without class assignment. nxq
#   alpha: Vector of hyperparameters for pi. Kx1
#   eta: Vector of hyperparameters for theta. dx1
#   mu0: List of vectors of mean hyperparameters for xi. List of K qx1 vectors
#   Sig0: List of matrices of variance hyperparameters for xi. List of K qxq matrices
# Outputs: `MCMC_out` list containing the following items:
#   pi_MCMC: Matrix of posterior samples for pi. (n_iter)xK
#   theta_MCMC: Array of posterior samples for theta. (n_iter)xpxKxd
#   xi_MCMC: Array of posterior samples for xi. (n_iter)xKxq
#   c_all_MCMC: Matrix of posterior samples for c_all. (n_iter)xn
#   z_all_MCMC: Matrix of posterior samples for z_all. (n_iter)xn
#   loglik_MCMC: Vector of posterior samples for log-likelihood. (n_iter)x1
run_MCMC_Rcpp <- function(OLCA_params, probit_params, n_runs, burn, thin, K, p, d, n, 
                          q, w_all, x_mat, y_all, V, alpha, eta, mu0, Sig0) {
  n_storage <- ceiling(n_runs / thin)  # Number of MCMC iterations to store
  
  # Initialize variables
  pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
  theta_MCMC <- array(NA, dim = c(n_storage, p, K, d))
  xi_MCMC <- array(NA, dim = c(n_storage, K, q))
  c_all_MCMC <- z_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
  loglik_MCMC <- numeric(n_storage)
  loglik <- numeric(n)  # Individual log-likelihood
  lin_pred <- numeric(n)                              # Linear predictor, V*xi
  
  
  # Initialized values
  pi <- OLCA_params$pi
  theta <- OLCA_params$theta
  c_all <- as.double(OLCA_params$c_all)  # allows updating by reference in rcpparmadillo
  xi <- probit_params$xi
  z_all <- probit_params$z_all
  
  # Update parameters and variables
  for (m in 1:n_runs) {
    update_pi(pi = pi, w_all = w_all, c_all = c_all, K = K, alpha = alpha)
    # print(round(pi, 3))
    update_c(c_all = c_all, n = n, K = K, p = p, theta = theta, 
             x_mat = x_mat, pi = pi, z_all = z_all, V = V, xi = xi, 
             y_all = y_all)
    # print(prop.table(table(c_all)))
    update_theta(theta = theta, p = p, K = K, d = d, eta = eta, 
                 w_all = w_all, c_all = c_all, x_mat = x_mat)
    # print(theta[1:5,1:10,1])
    # xi <- update_xi_R(n = n, w_all = w_all, V = V, K = K, c_all = c_all,
    #                   Sig0 = Sig0, mu0 = mu0, z_all = z_all, xi = xi)
    # z_all <- update_z_R(lin_pred = lin_pred, n = n, V = V, xi = xi,
    #                     c_all = c_all, z_all = z_all, y_all = y_all)
    xi <- update_xi(xi = xi, n = n, K = K, w_all = w_all, c_all = c_all,
                  z_all = z_all, V = V, mu0 = mu0, Sig0 = Sig0)
    # print(xi[1:3,])
    z_all <- update_z(z_all = z_all, n = n, V = V, xi = xi, c_all = c_all,
                    y_all = y_all)
    # print(z_all[1:10])
    update_loglik(loglik = loglik, n = n, p = p, c_all = c_all, 
                  theta = theta, x_mat = x_mat, pi = pi, 
                  z_all = z_all, V = V, xi = xi, y_all = y_all)
    # print(loglik[1:10])
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
  
  MCMC_out <- list(pi_MCMC = pi_MCMC, theta_MCMC = theta_MCMC, xi_MCMC = xi_MCMC,
                   c_all_MCMC = c_all_MCMC, z_all_MCMC = z_all_MCMC, 
                   loglik_MCMC = loglik_MCMC)
  return(MCMC_out)
}

# `get_mode` obtains the modal value given an input vector
# Inputs:
#   v: Input vector
# Outputs most common value found in input vector `v`
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# 'post_process' conducts post-processing to cluster individuals into a 
# reduced number of classes and reorder posterior parameter samples according to
# the reduced number of classes
# Inputs:
#   MCMC_out: Output from 'run_MCMC' containing 'pi_MCMC', 'theta_MCMC', 'xi_MCMC', 'c_all_MCMC', 'z_all_MCMC', and 'loglik_MCMC'
#   p: Number of exposure items
#   d: Number of exposure categories
#   q: Number of regression covariates excluding class assignment
# Outputs: `post_MCMC_out` list containing the following items:
#   K_med: Median, across iterations, of number of classes with >= 5% of individuals
#   pi: Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)
#   theta: Array of reduced and relabeled posterior samples for theta. (n_iter)xpx(K_med)xd
#   xi: Array of reduced and relabeled posterior samples for xi. (n_iter)x(K_med)xq
#   dendrogram: Hierarchical clustering dendrogram used for relabeling
post_process <- function(MCMC_out, p, d, q) {
  # Get median number of classes with >= 5% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
  
  # Cluster individuals into reduced number of classes using agglomerative clustering
  # Calculate pairwise distance matrix using Hamming distance: proportion of 
  # iterations where two individuals have differing class assignments
  distMat <- hamming.distance(t(MCMC_out$c_all_MCMC))
  dendrogram <- hclust(as.dist(distMat), method = "complete") # Hierarchical clustering dendrogram
  red_c_all <- cutree(dendrogram, k = K_med)                  # Group individuals into K_med classes
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

# `analyze_results` obtains posterior parameter samples and estimates
# Inputs:
#   MCMC_out: Output from 'run_MCMC' containing 'pi_MCMC', 'theta_MCMC', 'xi_MCMC', 
# 'c_all_MCMC', 'z_all_MCMC', 'loglik_MCMC'
#   post_MCMC_out: output from 'post_process' containing 'K_med', 'pi', 'theta', 'xi'
#   n: Number of individuals
#   p: Number of exposure items
#   V: Regression design matrix without class assignment. nxq
#   y_all: Vector of outcomes. nx1
#   x_mat: Categorical exposure matrix. nxp
# Outputs: `analysis` list containing the following items:
#   K_red: Number of unique classes
#   pi_red: Matrix of final posterior samples for pi. Mx(K_red)
#   theta_red: Array of final posterior samples for theta. Mxpx(K_red)xd
#   xi_red: Array of final posterior samples for xi. Mx(K_red)xq
#   pi_med: Vector of posterior median estimates for pi. (K_red)x1
#   theta_med: Array of posterior median estimates for theta. px(K_red)xd
#   xi_med: Matrix of posterior median estimates for xi. (K_red)xq
#   Phi_med: Vector of final individual outcome probabilities. nx1
#   c_all: Vector of final individual class assignments. nx1
#   pred_class_probs: Matrix of individual posterior class probabilities. nx(K_red)
#   loglik_med: Vector of final indiviudal log-likehoods. nx1
analyze_results <- function(MCMC_out, post_MCMC_out, n, p, V, y_all, x_mat) {
  
  #============== Identify unique classes using modal exposure categories ======
  # Posterior median estimate for theta across iterations
  theta_med_temp <- apply(post_MCMC_out$theta, c(2, 3, 4), median)
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
  theta_red <- aaply(theta_red, c(1, 2, 3), function(x) x / sum(x)) # Re-normalize
  xi_red <- post_MCMC_out$xi[, unique_classes, , drop = FALSE]
  
  #============== Posterior median estimates ===================================
  pi_med <- apply(pi_red, 2, median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), median, na.rm = TRUE)
  theta_med <- aaply(theta_med, c(1, 2), function(x) x / sum(x))  # Re-normalize
  xi_med <- apply(xi_red, c(2, 3), median, na.rm = TRUE)
  
  #============== Update c using unique classes and posterior estimates ========
  z_all <- MCMC_out$z_all_MCMC[M, ]
  c_all <- MCMC_out$c_all_MCMC[M, ]                      # Final class assignments
  pred_class_probs <- matrix(NA, nrow = n, ncol = K_red) # Posterior class membership probabilities
  log_cond_c <- matrix(NA, nrow = n, ncol = K_red)       # Individual log-likelihood for each class
  Phi_med_all_c <- pnorm(V %*% t(xi_med))  # Outcome probabilities for all classes
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
      log_probit_part <- log(dnorm(z_all[i], mean = V[i, ] %*% xi_med[k, ])) 
      if (log_probit_part == -Inf) {
        log_probit_part <- log(1e-16)
      }
      log_probit_comp_k <- log_probit_part + 
        log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
      # Individual log-likelihood for class k
      log_cond_c[i, k] <- log(pi_med[k]) + log_theta_comp_k + log_probit_comp_k
    }
    # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs[i, ] <- exp(log_cond_c[i, ] - logSumExp(log_cond_c[i, ]))
    # Update class assignment using the posterior probabilities
    c_all[i] <- rcat(n = 1, p = pred_class_probs[i, ])
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
      log(dnorm(z_all[i], mean = V[i, ] %*% xi_med[c_i, ])) + 
      log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
  }
  
  analysis <- list(K_red = K_red, pi_red = pi_red, theta_red = theta_red, 
                   xi_red = xi_red, pi_med = pi_med, theta_med = theta_med, 
                   xi_med = xi_med, Phi_med = Phi_med, c_all = c_all, 
                   pred_class_probs = pred_class_probs, loglik_med = loglik_med)
  return(analysis)
}

# `unconstrain` is a helper function that converts each row of an input array of 
# MCMC parameter output from constrained space to unconstrained space in Stan
# Inputs:
#   i: Row index
#   K: Number of classes
#   stan_model: Stan model
#   pi: MCMC matrix output for pi; MxK
#   theta: MCMC array output for theta; MxpxKxd
#   xi: MCMC matrix output for xi; MxKxS
# Outputs vector of unconstrained parameters
unconstrain <- function(i, K, stan_model, pi, theta, xi) {
  # Be careful with dimension of xi when latent class is only covariate, as R
  # will automatically drop dimensions of size 1
  u_pars <- unconstrain_pars(stan_model, 
                             list("pi" = pi[i,], "theta" = theta[i,,,], 
                                  "xi" = as.matrix(xi[i,,])))
  return(u_pars)
}

# `DEadj` is a helper function to apply the matrix rotation
# Inputs:
#   par: Unadjusted parameter estimates
#   par_hat: Unadjusted median parameter estimates
#   R2R1: Adjustment matrix
# Outputs adjusted parameter estimates
DEadj <- function(par, par_hat, R2R1) {
  par_adj <- (par - par_hat) %*% R2R1 + par_hat
  par_adj <- as.vector(par_adj)
  return(par_adj)
}

# `grad_par` is a helper function nested in 'withReplicates()' to obtain the 
# gradient with the replicate weights
# Inputs:
#   pwts: Replicate weights from 'svyrepdesign' object
#   svydata: Data frame containing all variables from 'svyrepdesign' object
#   stan_mod: Stan model object
#   stan_data: Stan data input
#   par_stan: Parameters with respect to which gradient should be computed
#   u_pars: Unconstrained parameters estimates for evaluating gradient
# Outputs `gradpar` gradient evaluated at `u_pars` using replicate weights
grad_par <- function(pwts, svydata, stan_mod, stan_data, par_stan, u_pars) {
  stan_data$weights <- pwts
  out_stan <- sampling(object = stan_mod, data = stan_data, pars = par_stan,
                       chains = 0, iter = 0, refresh = 0)
  gradpar <- grad_log_prob(out_stan, u_pars)
  return(gradpar)
}

# `var_adj` applies applies the post-processing variance adjustment
# Inputs: 
#   mod_stan: Stan model
#   analysis: Output from `analyze_results` containing K_red, pi_red, theta_red, 
# xi_red, pi_med, theta_med, xi_med, Phi_med, c_all, pred_class_probs, loglik_med
#   K: Number of classes
#   p: Number of exposure items
#   d: Number of exposure categories
#   n: Number of individuals
#   q: Number of regression covariates excluding class assignment
#   x_mat: Categorical exposure matrix. nxp
#   y_all: Vector of outcomes. nx1
#   V: Regression design matrix without class assignment. nxq
#   w_all: Weights normalized to sum to n. nx1
#   s_all: Vector of stratifying variable for individuals. nx1
#   clus_id_all: Vector of cluster indicators. nx1 
# Outputs: `analysis_adj` list containing the following items:
#   pi_red_adj: Matrix of adjusted posterior samples for pi. Mx(K_red)
#   theta_red_adj: Array of adjusted posterior samples for theta. Mxpx(K_red)xd
#   xi_red_adj: Array of adjusted posterior samples for xi. Mx(K_red)xq
#   pi_med_adj: Vector of adjusted posterior median estimates for pi. (K_red)x1
#   theta_med_adj: Array of adjusted posterior median estimates for theta. px(K_red)xd
#   xi_med_adj: Matrix of adjusted posterior median estimates for xi. (K_red)xq
#   Phi_med: Vector of adjusted individual outcome probabilities. nx1
#   c_all: Vector of final individual class assignments from `analyze_results`. nx1
#   pred_class_probs: Matrix of individual posterior class probabilities from `analyze_results`. nx(K_red)
#   loglik_med: Vector of final individual log-likelihoods from `analyze_results`. nx1
var_adjust <- function(mod_stan, analysis, K, p, d, n, q, x_mat, y_all, V, w_all, 
                       s_all, clus_id_all) {
  #=============== Run Stan model ==============================================
  # Define data for Stan model
  alpha <- rep(1, K) / K            # Hyperparameter for prior for pi
  eta <- rep(1, d)                  # Hyperparameter for prior for theta
  mu0 <- Sig0 <- vector("list", K)  # Hyperparameters for xi
  for (k in 1:K) {
    mu0[[k]] <- rnorm(n = q)
    Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
  }
  data_stan <- list(K = K, p = p, d = d, n = n, q = q, X = x_mat, y = y_all, 
                    V = V, weights = w_all, alpha = alpha, eta = eta, mu0 = mu0, 
                    Sig0 = Sig0)
  
  # Stan parameters of interest
  par_stan <- c('pi', 'theta', 'xi')  # subset of parameters interested in
  
  # Run Stan model
  # Stan will pass warnings from calling 0 chains, but will still create an 
  # out_stan object for the 'grad_log_prob()' method
  out_stan <- sampling(object = mod_stan, data = data_stan, pars = par_stan,
                       chains = 0, iter = 0, refresh = 0)
  
  #=============== Convert to unconstrained parameters =========================
  # Convert params from constrained space to unconstrained space
  unc_par_hat <- unconstrain_pars(out_stan, list("pi" = analysis$pi_med,
                                                 "theta" = analysis$theta_med,
                                                 "xi" = analysis$xi_med))
  # Get posterior MCMC samples in unconstrained space for all parameters
  M <- dim(analysis$pi_red)[1]
  unc_par_samps <- lapply(1:M, unconstrain, stan_model = out_stan, K = K, 
                          pi = analysis$pi_red, theta = analysis$theta_red, 
                          xi = analysis$xi_red)
  unc_par_samps <- matrix(unlist(unc_par_samps), byrow = TRUE, nrow = M)
  
  #=============== Post-processing adjustment in unconstrained space ===========
  # Estimate Hessian
  H_hat <- -1*optimHess(unc_par_hat, gr = function(x){grad_log_prob(out_stan, x)})
  
  # Create survey design
  if (!is.null(s_all)) {  # Include stratifying variable
    svy_data <- data.frame(s = s_all, 
                           x = x_mat,
                           y = y_all, 
                           wts = w_all,
                           clus = clus_id_all)
    svydes <- svydesign(ids = ~clus, strata = ~s, weights = ~wts, data = svy_data)
  } else {  # No stratifying variable
    svy_data <- data.frame(x = x_mat,
                           y = y_all, 
                           wts = w_all,
                           clus = clus_id_all)
    svydes <- svydesign(ids = ~clus, weights = ~wts, data = svy_data)
  }
  
  # Create svrepdesign
  svyrep <- as.svrepdesign(design = svydes, type = "mrbbootstrap", 
                           replicates = 100)
  
  rep_temp <- withReplicates(design = svyrep, theta = grad_par, 
                             stan_mod = mod_stan, stan_data = data_stan, 
                             par_stan = par_stan, u_pars = unc_par_hat)
  J_hat <- vcov(rep_temp)
  
  # Compute adjustment
  H_inv <- solve(H_hat)
  V1 <- H_inv %*% J_hat %*% H_inv
  
  # Check for issues with negative diagonals
  if (min(diag(V1)) < 0) {
    print("V1 has negative variances")
  }
  if (min(diag(H_inv)) < 0) {
    print("H_inv has negative variances")
  }
  # If matrices are not p.d. due to rounding issues, convert to nearest p.d. 
  # matrix using method proposed in Higham (2002)
  if (min(Re(eigen(V1)$values)) < 0) { 
    V1_pd <- nearPD(V1)
    R1 <- chol(V1_pd$mat)
    print(paste0("V1: absolute eigenvalue difference to nearest p.d. matrix: ", 
                 sum(abs(eigen(V1)$values - eigen(V1_pd$mat)$values))))
  } else {
    R1 <- chol(V1)
  }
  if (min(Re(eigen(H_inv)$values)) < 0) {
    H_inv_pd <- nearPD(H_inv)
    R2_inv <- chol(H_inv_pd$mat)
    print(paste0("H_inv: absolute eigenvalue difference to nearest p.d. matrix: ", 
                 sum(abs(eigen(H_inv)$values - eigen(H_inv_pd$mat)$values))))
  } else {
    R2_inv <- chol(H_inv)
  }
  # Obtain the variance adjustment matrix
  R2 <- solve(R2_inv)
  R2R1 <- R2 %*% R1
  
  # Apply variance adjustment to parameters
  par_adj <- apply(unc_par_samps, 1, DEadj, par_hat = unc_par_hat, R2R1 = R2R1, 
                   simplify = FALSE)
  par_adj <- matrix(unlist(par_adj), byrow=TRUE, nrow=M)
  
  #=============== Convert adjusted to constrained space =======================
  # Constrained adjusted parameters for all MCMC samples
  pi_red_adj <- matrix(NA, nrow=M, ncol=K)
  theta_red_adj <- array(NA, dim=c(M, p, K, d))
  xi_red_adj <- array(NA, dim=c(M, K, q))
  for (i in 1:M) {
    constr_pars <- constrain_pars(out_stan, par_adj[i,])
    pi_red_adj[i, ] <- constr_pars$pi
    theta_red_adj[i,,,] <- constr_pars$theta
    xi_red_adj[i,,] <- constr_pars$xi
  }
  
  #=============== Output adjusted parameters ==================================
  pi_med_adj <- apply(pi_red_adj, 2, median)
  theta_med_adj <- apply(theta_red_adj, c(2,3,4), median)
  xi_med_adj <- apply(xi_red_adj, c(2,3), median)
  
  # Update Phi_med using adjusted xi estimate
  c_all <- analysis$c_all
  Phi_med_all_c <- pnorm(V %*% t(xi_med_adj))  # Outcome probabilities for all classes
  Phi_med_adj <- numeric(n)                    # Initialize individual outcome probabilities
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    # Calculate outcome probabilities P(Y=1|-) using updated class assignment
    Phi_med_adj[i] <- Phi_med_all_c[i, c_all[i]] 
  }
  
  analysis_adj <- list(pi_red_adj = pi_red_adj, theta_red_adj = theta_red_adj, 
                       xi_red_adj = xi_red_adj, pi_med_adj = pi_med_adj, 
                       theta_med_adj = theta_med_adj, xi_med_adj = xi_med_adj, 
                       Phi_med_adj = Phi_med_adj, c_all = c_all,
                       pred_class_probs = analysis$pred_class_probs,
                       log_lik_med = analysis$loglik_med)
  return(analysis_adj)
}

# `run_MCMC_Rcpp_WOLCA` runs the Gibbs sampler MCMC algorithm to obtain posterior 
# samples for the two-step unsupervised WOLCA model
# Inputs:
#   OLCA_params: output from 'init_OLCA' containing 'pi', 'c_all', 'theta'
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
#   K: Number of classes
#   p: Number of exposure items
#   d: Number of exposure categories
#   n: Number of individuals
#   w_all: Weights normalized to sum to n. nx1
#   x_mat: Categorical exposure matrix. nxp
#   alpha: Vector of hyperparameters for pi. Kx1
#   eta: Vector of hyperparameters for theta. dx1
# Outputs: `MCMC_out` list containing the following items:
#   pi_MCMC: Matrix of posterior samples for pi. (n_iter)xK
#   theta_MCMC: Array of posterior samples for theta. (n_iter)xpxKxd
#   c_all_MCMC: Matrix of posterior samples for c_all. (n_iter)xn
run_MCMC_Rcpp_WOLCA <- function(OLCA_params, n_runs, burn, thin, K, p, d, n, 
                                w_all, x_mat, alpha, eta) {
  n_storage <- ceiling(n_runs / thin)  # Number of MCMC iterations to store
  
  # Initialize variables
  pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
  theta_MCMC <- array(NA, dim = c(n_storage, p, K, d))
  c_all_MCMC <- z_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
  
  # Initialized values
  pi <- OLCA_params$pi
  theta <- OLCA_params$theta
  c_all <- OLCA_params$c_all
  
  # Update parameters and variables
  for (m in 1:n_runs) {
    update_pi(pi = pi, w_all = w_all, c_all = c_all, K = K, alpha = alpha)
    # print(pi[1:10])
    update_c_WOLCA(c_all = c_all, n = n, K = K, p = p, theta = theta, 
                            x_mat = x_mat, pi = pi)
    # print(c_all[1:10])
    update_theta(theta = theta, p = p, K = K, d = d, eta = eta, 
                         w_all = w_all, c_all = c_all, x_mat = x_mat)
    # print(theta[1:5,1:10,1])
    
    #============== Store posterior values based on thinning  ==================
    if (m %% thin == 0) {
      m_thin <- m / thin
      pi_MCMC[m_thin, ] <- pi
      theta_MCMC[m_thin, , , ] <- theta
      c_all_MCMC[m_thin, ] <- c_all
    }
    
    #============== Relabel classes every 10 iterations to encourage mixing ====
    if (m %% 10 == 0) {
      new_order <- permute(1:K)      # New permuted labels
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

# 'post_process_WOLCA' conducts post-processing to cluster individuals into a 
# reduced number of classes and reorder posterior parameter samples according to
# the reduced number of classes
# Inputs:
#   MCMC_out: Output from 'run_MCMC' containing 'pi_MCMC', 'theta_MCMC', 'xi_MCMC', 'c_all_MCMC', 'z_all_MCMC', and 'loglik_MCMC'
#   p: Number of exposure items
#   d: Number of exposure categories
# Outputs: `post_MCMC_out` list containing the following items:
#   K_med: Median, across iterations, of number of classes with >= 5% of individuals
#   pi: Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)
#   theta: Array of reduced and relabeled posterior samples for theta. (n_iter)xpx(K_med)xd
post_process_WOLCA <- function(MCMC_out, p, d) {
  # Get median number of classes with >= 5% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
  
  # Cluster individuals into reduced number of classes using agglomerative clustering
  # Calculate pairwise distance matrix using Hamming distance: proportion of 
  # iterations where two individuals have differing class assignments
  distMat <- hamming.distance(t(MCMC_out$c_all_MCMC))
  dendrogram <- hclust(as.dist(distMat), method = "complete") # Hierarchical clustering dendrogram
  red_c_all <- cutree(dendrogram, k = K_med)                  # Group individuals into K_med classes
  # For each iteration, relabel new classes using the most common old class assignment
  relabel_red_classes <- matrix(NA, nrow = M, ncol = K_med)   # Apply new classes to each iteration
  for (k in 1:K_med) {
    relabel_red_classes[, k] <- apply(as.matrix(MCMC_out$c_all_MCMC[, red_c_all == k]), 
                                      1, get_mode)
  }
  
  # Reduce and reorder parameter estimates using new classes
  pi <- matrix(NA, nrow = M, ncol = K_med)
  theta <- array(NA, dim = c(M, p, K_med, d))
  for (m in 1:M) {
    iter_order <- relabel_red_classes[m, ]
    pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
    pi[m, ] <- pi_temp / sum(pi_temp)
    theta[m, , , ] <- MCMC_out$theta_MCMC[m, , iter_order, ]
  }
  
  post_MCMC_out <- list(K_med = K_med, pi = pi, theta = theta)
  
  return(post_MCMC_out)  
  # plot(dendrogram, labels = FALSE)
  # rect.hclust(dendrogram, k = K_med)
}

# `analyze_results_WOLCA` obtains posterior parameter samples and estimates
# Inputs:
#   MCMC_out: Output from 'run_MCMC_Rcpp_WOLCA' containing 'pi_MCMC', 
# 'theta_MCMC', 'c_all_MCMC'
#   post_MCMC_out: output from 'post_process_WOLCA' containing 'K_med', 'pi', 'theta'
#   n: Number of individuals
#   p: Number of exposure items
#   x_mat: Categorical exposure matrix. nxp
# Outputs: `analysis` list containing the following items:
#   K_red: Number of unique classes
#   pi_red: Matrix of final posterior samples for pi. Mx(K_red)
#   theta_red: Array of final posterior samples for theta. Mxpx(K_red)xd
#   pi_med: Vector of posterior median estimates for pi. (K_red)x1
#   theta_med: Array of posterior median estimates for theta. px(K_red)xd
#   c_all: Vector of final individual class assignments. nx1
#   pred_class_probs: Matrix of individual posterior class probabilities. nx(K_red)
analyze_results_WOLCA <- function(MCMC_out, post_MCMC_out, n, p, x_mat) {
  
  #============== Identify unique classes using modal exposure categories ======
  # Posterior median estimate for theta across iterations
  theta_med_temp <- apply(post_MCMC_out$theta, c(2, 3, 4), median)
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
  theta_red <- aaply(theta_red, c(1, 2, 3), function(x) x / sum(x)) # Re-normalize
  
  #============== Posterior median estimates ===================================
  pi_med <- apply(pi_red, 2, median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), median, na.rm = TRUE)
  theta_med <- aaply(theta_med, c(1, 2), function(x) x / sum(x))  # Re-normalize
  
  #============== Update c using unique classes and posterior estimates ========
  c_all <- MCMC_out$c_all_MCMC[M, ]                      # Final class assignments
  pred_class_probs <- matrix(NA, nrow = n, ncol = K_red) # Posterior class membership probabilities
  log_cond_c <- matrix(NA, nrow = n, ncol = K_red)       # Individual log-likelihood for each class
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    for (k in 1:K_red) {
      # Calculate theta component of individual log-likelihood assuming class k
      log_theta_comp_k <- 0
      for (j in 1:p) {
        log_theta_comp_k <- log_theta_comp_k + log(theta_med[j, k, x_mat[i, j]])
      }
      # Individual log-likelihood for class k
      log_cond_c[i, k] <- log(pi_med[k]) + log_theta_comp_k
    }
    # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs[i, ] <- exp(log_cond_c[i, ] - logSumExp(log_cond_c[i, ]))
    # Update class assignment using the posterior probabilities
    c_all[i] <- rcat(n = 1, p = pred_class_probs[i, ])
  }
  
  analysis <- list(K_red = K_red, pi_red = pi_red, theta_red = theta_red, 
                   pi_med = pi_med, theta_med = theta_med, c_all = c_all, 
                   pred_class_probs = pred_class_probs)
  
  return(analysis)
}

# `manual_CI` manually calculates a Wald confidence interval using a t-dist
# with df from the survey design. Used for WOLCA in the situation where svyglm 
# produces negative residual df, calculated as design df plus one, minus the 
# number of parameters estimated. Best if no cluster-level covariates in the 
# regression model
# Inputs:
#   model_object: svyglm model fit object
#   svy_df: survey design df
#   ci: confidence interval level
# Outputs: dataframe of confidence interval for all coefficients
manual_CI <- function(model_object, svy_df, ci = 0.95){
  a <- coef(summary(model_object))
  mult <- qt((1 + ci) / 2, df = svy_df)
  restab <- with(as.data.frame(a),
                 cbind(est = Estimate,
                       lwr =  Estimate - mult*`Std. Error`,
                       upr = Estimate + mult*`Std. Error`))
  rownames(restab) <- rownames(a)
  return(data.frame(restab))
}


# `convert_ref_to_comb` converts a combination of factor variable and 
# reference cell coding to purely reference cell coding
# Inputs:
#   beta_ref: Matrix of probit coefficients in reference cell coding; (K*q)x1
# Outputs:
#   beta_comb: Matrix of probit coefficients in combination coding; Kxq
convert_ref_to_comb <- function(beta_ref) {
  beta_comb <- beta_ref
  for (i in 2:nrow(beta_ref)) {
    beta_comb[i, ] <- beta_ref[i, ] - beta_ref[1, ]
  }
  return(beta_comb)
}
