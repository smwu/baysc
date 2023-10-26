#' Helper function that converts each row of an input array of MCMC parameter 
#' output from constrained space to unconstrained space in Stan
#' 
#' @param i Row index
#' @param K Number of classes
#' @param stan_model Stan model
#' @param pi MCMC matrix output for pi; MxK
#' @param theta MCMC array output for theta; MxpxKxd
#' @param xi: MCMC matrix output for xi; MxKxS
#' 
#' @return Outputs vector of unconstrained parameters
#' @noRd
unconstrain <- function(i, K, stan_model, pi, theta, xi) {
  # Be careful with dimension of xi when latent class is only covariate, as R
  # will automatically drop dimensions of size 1
  u_pars <- unconstrain_pars(stan_model, 
                             list("pi" = pi[i,], "theta" = theta[i,,,], 
                                  "xi" = as.matrix(xi[i,,])))
  return(u_pars)
}

#' Helper function to apply the matrix rotation
#' 
#' @param par Unadjusted parameter estimates
#' @param par_hat Unadjusted median parameter estimates
#' @param R2R1 Adjustment matrix
#' 
#' @return Outputs adjusted parameter estimates
#' @noRd
DEadj <- function(par, par_hat, R2R1) {
  par_adj <- (par - par_hat) %*% R2R1 + par_hat
  par_adj <- as.vector(par_adj)
  return(par_adj)
}

#' Helper function nested in `withReplicates()` to obtain the gradient with the 
#' replicate weights
#' 
#' @param pwts Replicate weights from `svyrepdesign` object
#' @param svydata Data frame containing all variables from `svyrepdesign` object
#' @param stan_mod Stan model object
#' @param stan_data Stan data input
#' @param par_stan Parameters with respect to which gradient should be computed
#' @param u_pars Unconstrained parameters estimates for evaluating gradient
#' 
#' @return Outputs `gradpar` gradient evaluated at `u_pars` using replicate weights
#' @noRd
grad_par <- function(pwts, svydata, stan_mod, stan_data, par_stan, u_pars) {
  stan_data$weights <- pwts
  out_stan <- sampling(object = stan_mod, data = stan_data, pars = par_stan,
                       chains = 0, iter = 0, refresh = 0)
  gradpar <- grad_log_prob(out_stan, u_pars)
  return(gradpar)
}

#' Post-processing variance adjustment
#' 
#' @description 
#' `var_adj` applies applies the post-processing variance adjustment
#' 
#' @inheritParams run_MCMC_Rcpp
#' @param mod_stan Stan model
#' @param estimates Output from `get_estimates()` containing `K_red`, `pi_red`, 
#' `theta_red`, `xi_red`, `pi_med`, `theta_med`, `xi_med`, `Phi_med`, `c_all`, 
#' `pred_class_probs`, `loglik_med`
#' @param stratum_id Vector of stratifying variable for individuals. nx1
#' @param cluster_id Vector of cluster indicators. nx1 
#' 
#' @return 
#' Returns list `estimates_adj` containing:
#' 
#' \describe{
#'   \item{\code{pi_red_adj}}{Matrix of adjusted posterior samples for pi. Mx(K_red)}
#'   \item{\code{theta_red_adj}}{Array of adjusted posterior samples for theta. Mxpx(K_red)xd}
#'   \item{\code{xi_red_adj}}{Array of adjusted posterior samples for xi. Mx(K_red)xq}
#'   \item{\code{pi_med_adj}}{Vector of adjusted posterior median estimates for pi. (K_red)x1}
#'   \item{\code{theta_med_adj}}{Array of adjusted posterior median estimates for theta. px(K_red)xd}
#'   \item{\code{xi_med_adj}}{Matrix of adjusted posterior median estimates for xi. (K_red)xq}
#'   \item{\code{Phi_med_adj}}{Vector of adjusted individual outcome probabilities. nx1}
#'   \item{\code{c_all}}{Vector of final individual class assignments from `get_estimates()`. nx1}
#'   \item{\code{pred_class_probs}}{Matrix of individual posterior class probabilities from `get_estimates()`. nx(K_red)}
#'   \item{\code{loglik_med}}{Vector of final indiviudal log-likehoods from `get_estimates()`. nx1} 
#' }
#'
#' @export
#'
#' #examples
#' 
var_adjust <- function(mod_stan, estimates, K, p, d, n, q, x_mat, y_all, V, w_all, 
                       stratum_id, cluster_id) {
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
  unc_par_hat <- unconstrain_pars(out_stan, list("pi" = estimates$pi_med,
                                                 "theta" = estimates$theta_med,
                                                 "xi" = estimates$xi_med))
  # Get posterior MCMC samples in unconstrained space for all parameters
  M <- dim(estimates$pi_red)[1]
  unc_par_samps <- lapply(1:M, unconstrain, stan_model = out_stan, K = K, 
                          pi = estimates$pi_red, theta = estimates$theta_red, 
                          xi = estimates$xi_red)
  unc_par_samps <- matrix(unlist(unc_par_samps), byrow = TRUE, nrow = M)
  
  #=============== Post-processing adjustment in unconstrained space ===========
  # Estimate Hessian
  H_hat <- -1*optimHess(unc_par_hat, gr = function(x){grad_log_prob(out_stan, x)})
  
  # Create survey design
  if (!is.null(stratum_id)) {  # Include stratifying variable
    svy_data <- data.frame(s = stratum_id, 
                           x = x_mat,
                           y = y_all, 
                           wts = w_all,
                           clus = cluster_id)
    svydes <- svydesign(ids = ~clus, strata = ~s, weights = ~wts, data = svy_data)
  } else {  # No stratifying variable
    svy_data <- data.frame(x = x_mat,
                           y = y_all, 
                           wts = w_all,
                           clus = cluster_id)
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
  c_all <- estimates$c_all
  Phi_med_all_c <- pnorm(V %*% t(xi_med_adj))  # Outcome probabilities for all classes
  Phi_med_adj <- numeric(n)                    # Initialize individual outcome probabilities
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    # Calculate outcome probabilities P(Y=1|-) using updated class assignment
    Phi_med_adj[i] <- Phi_med_all_c[i, c_all[i]] 
  }
  
  estimates_adj <- list(pi_red_adj = pi_red_adj, theta_red_adj = theta_red_adj, 
                       xi_red_adj = xi_red_adj, pi_med_adj = pi_med_adj, 
                       theta_med_adj = theta_med_adj, xi_med_adj = xi_med_adj, 
                       Phi_med_adj = Phi_med_adj, c_all = c_all,
                       pred_class_probs = estimates$pred_class_probs,
                       log_lik_med = estimates$loglik_med)
  return(estimates_adj)
}