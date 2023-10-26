#===================================================
## Helper functions for WSOLCA, SOLCA, and WOLCA
## Programmer: SM Wu   
## Data: Simulations and application   
#===================================================


# `get_mode` obtains the modal value given an input vector
# Input: 
#   v: Input vvector
# Outputs most common value found in input vector `v`
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
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


