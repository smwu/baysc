#' Unconstrain parameters
#'
#' `unconstrain` is a helper function that converts each row of an input array 
#' of MCMC parameter output from constrained space to unconstrained space in Stan
#' 
#' @param i Row index
#' @param K Number of classes
#' @param stan_model Stan model
#' @param pi MCMC matrix output for pi; MxK
#' @param theta MCMC array output for theta; MxJxKxR
#' @param xi: MCMC matrix output for xi; MxKxS
#' 
#' @return Outputs vector of unconstrained parameters
#' @importFrom rstan unconstrain_pars
#' @keywords internal
#' @export
unconstrain <- function(i, K, stan_model, pi, theta, xi) {
  # Be careful with dimension of xi when latent class is only covariate, as R
  # will automatically drop dimensions of size 1
  u_pars <- rstan::unconstrain_pars(stan_model, 
                             list("pi" = pi[i,], "theta" = theta[i,,,], 
                                  "xi" = as.matrix(xi[i,,])))
  return(u_pars)
}

#' Adjust estimates
#'
#' @description
#' `DEadj` is a helper function to apply the matrix rotation
#' 
#' @param par Unadjusted parameter estimates
#' @param par_hat Unadjusted median parameter estimates
#' @param R2R1 Adjustment matrix
#' 
#' @return Outputs adjusted parameter estimates
#' @keywords internal
#' @export
DEadj <- function(par, par_hat, R2R1) {
  par_adj <- (par - par_hat) %*% R2R1 + par_hat
  par_adj <- as.vector(par_adj)
  return(par_adj)
}

#' Get gradient of log posterior from the unconstrained parameter space
#' 
#' @description
#' `grad_par` is a helper function nested in `withReplicates()` to obtain the 
#' gradient with the replicate weights
#' 
#' @param pwts Replicate weights from `svyrepdesign` object
#' @param svydata Data frame containing all variables from `svyrepdesign` object
#' @param stan_mod Stan model object
#' @param stan_data Stan data input
#' @param par_stan Parameters with respect to which gradient should be computed
#' @param u_pars Unconstrained parameters estimates for evaluating gradient
#' 
#' @importFrom rstan sampling grad_log_prob
#' @return Outputs `gradpar` gradient evaluated at `u_pars` using replicate weights
#' @keywords internal
#' @export
grad_par <- function(pwts, svydata, stan_mod, stan_data, par_stan, u_pars) {
  stan_data$weights <- pwts
  out_stan <- rstan::sampling(object = stan_mod, data = stan_data, pars = par_stan,
                       chains = 0, iter = 0, refresh = 0)
  gradpar <- rstan::grad_log_prob(out_stan, u_pars)
  return(gradpar)
}


#' Post-processing variance adjustment
#' 
#' @description 
#' `swolca_var_adj` applies applies the post-processing variance adjustment after 
#' a call to [swolca()] to correct for underestimation of posterior intervals.
#' 
#' @param res An object of class `"swolca"`, resulting from a call to [swolca()],
#' containing the unadjusted estimates.
#' @param alpha Hyperparameter for prior for class membership probabilities 
#' \eqn{\pi}. Default is `NULL` and default values are used (see Details below). 
#' If specified, must be (`K_max`)x1. 
#' @param eta Hyperparameter for prior for item consumption level probabilities 
#' \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class \eqn{k}, assumed to be 
#' the same across classes. Default is `NULL` and default values are used (see 
#' Details below). If specified, must be JxR, where J is the number of exposure 
#' items and R is the maximum number of levels for the exposure.
#' @param mu0 Mean hyperparameters for regression coefficients \eqn{\xi_{k\cdot}} 
#' for each class \eqn{k}. Default is `NULL` and default values are used (see 
#' Details below). If specified, must be a list of `K_max` vectors of dimension 
#' Qx1, where Q is the number of regression covariates excluding latent class assignment.  
#' @param Sig0 Variance hyperparameters for regression coefficients 
#' \eqn{\xi_{k\cdot}} for each class \eqn{k}. Default is `NULL` and default 
#' values are used (see Details below). If specified, must be a list of `K_max` 
#' QxQ matrices, where Q is the number of regression covariates excluding latent 
#' class assignment. 
#' @param num_reps Number of bootstrap replicates to use for the variance 
#' adjustment estimate. Default is 100.
#' @param save_res Boolean specifying if results should be saved. Default = `TRUE`.
#' @param save_path String specifying directory and file name to save results, 
#' e.g., "~/Documents/run". Default is `NULL`. If this is the same as the file 
#' name specified in [swolca()], the unadjusted results are overwritten with the 
#' adjusted results.
#' @param adjust_seed Numeric seed for variance adjustment. Default is `NULL`.
#' 
#' @details
#' `var_adjust` applies a post-processing variance adjustment that rescales the
#' variance to obtain correct coverage of posterior intervals, adapted from 
#' Williams and Savitsky (2021). To obtain the rescaling, a sandwich-type 
#' variance is estimated. To estimate the Hessian that composes the "bread" of 
#' the sandwich, the mixture model is specified in Stan and the parameters are 
#' converted to the unconstrained space. Bootstrap replicates are used to 
#' estimate the covariance matrix that composes the "meat" of the sandwich. 
#' If there are any rounding issues, the resulting matrices are mapped to the 
#' nearest positive definite matrix. Next, the rescaling adjustment is derived
#' and applied to the parameters for all MCMC iterations. Finally, the adjusted
#' parameters are converted back to the constrained space and the posterior 
#' median parameter estimates are recomputed, now with proper variance estimation. 
#' 
#' To save results, set `save_res = TRUE` (default) and `save_path` to a string
#' that specifies both the location and the beginning of the file name 
#' (e.g., "~/Documents/run"). The file name will have "_swolca_results.RData" 
#' appended to it, overwriting the unadjusted results if the file names are the 
#' same.
#' 
#' If hyperparameters are left as `NULL` (default), the following default 
#' values are used. Let \eqn{K} refer to the final number of latent class 
#' obtained from running [swolca()], available at `res$estimates$K_red`.
#' For \eqn{\pi}, a Dirichlet prior with hyperparameter \eqn{\alpha = 1/K} for 
#' each component. For \eqn{\theta_{jk\cdot}}, a Dirichlet prior with 
#' hyperparameter  \eqn{\eta_j} equal to `rep(1, R_j)` where `R_j` is the number 
#' of levels for exposure item j. If `R_j < R`, the remaining levels have
#' hyperparameter set to 0.01. This is done independently for each exposure item j
#' and is assumed to be the same across latent classes. For \eqn{\xi_{k\cdot}}, a 
#' Multivariate Normal distribution with mean vector hyperparameter \eqn{\mu_0} 
#' drawn from a Normal(0,1) hyperprior for each component, and variance matrix 
#' hyperparameter \eqn{\Sigma_0} a diagonal matrix with diagonal components drawn 
#' from InvGamma(shape=3.5, scale=6.25) distributions. 
#' 
#' If the following warning message appears, please run the sampler for more 
#' iterations as there is instability in the parameter estimates, usually in the 
#' form of large positive or negative values for \eqn{\xi}): 
#' "Error in svrVar(thetas, scale, rscales, mse = design$mse, coef = full) : 
#' All replicates contained NAs".
#' 
#' When running the variance adjustment, the following warning messages may appear: 
#' "the number of chains is less than 1; sampling not done" and 
#' "In mrbweights(design$cluster, design$strata, design$fpc, ...) : Design is 
#' sampled with replacement: only first stage used". These messages do not pose 
#' an issue to the statistical validity of the methods and can be ignored. 
#' 
#' @return 
#' Returns an object `res` of class `"swolca"`, which includes all outputs from 
#' [swolca()] as well as a list `estimates_adjust` containing:
#' \describe{
#'   \item{\code{pi_red}}{Matrix of adjusted posterior samples for pi. Mx(K_red), 
#'   where M is the number of MCMC iterations after burn-in and thinning.}
#'   \item{\code{theta_red}}{Array of adjusted posterior samples for theta. MxJx(K_red)xR}
#'   \item{\code{xi_red}}{Array of adjusted posterior samples for xi. Mx(K_red)xQ}
#'   \item{\code{pi_med}}{Vector of adjusted posterior median estimates for pi. (K_red)x1}
#'   \item{\code{theta_med}}{Array of adjusted posterior median estimates for theta. px(K_red)xR}
#'   \item{\code{xi_med}}{Matrix of adjusted posterior median estimates for xi. (K_red)xQ}
#'   \item{\code{Phi_med}}{Vector of adjusted individual outcome probabilities. nx1}
#'   \item{\code{c_all}}{Vector of final individual class assignments from `swolca()`. nx1}
#'   \item{\code{pred_class_probs}}{Matrix of individual posterior class 
#'   probabilities from `swolca()`. nx(K_red)}
#'   \item{\code{loglik_med}}{Vector of final indiviudal log-likehoods from `swolca()`. nx1} 
#' }
#' The `runtime` output for `res` is also updated to include the runtime for the 
#' variance adjustment in addition to the runtime for the main `swolca()` model.
#' 
#' If `save_res = TRUE` (default), the updated `res` object is saved as 
#' `[save_path]_swolca_results.RData`, overwriting the unadjusted results if the 
#' file names are the same. 
#'
#' @seealso [swolca()]
#' @importFrom stats rnorm pnorm optimHess vcov median
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom rstan sampling unconstrain_pars grad_log_prob constrain_pars
#' @importFrom survey svydesign as.svrepdesign withReplicates
#' @importFrom Matrix nearPD
#' @export
#' 
#' @references 
#' Williams, M. R., & Savitsky, T. D. (2021). Uncertainty Estimation for 
#' Pseudo‐Bayesian Inference Under Complex Sampling. International Statistical 
#' Review, 89(1), 72-107.
#'
#' @examples
#' \dontrun{   
#' # Load simulated data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
#' stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
#' sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
#' n <- dim(x_mat)[1]                   # Number of individuals
#' 
#' # Probit model only includes latent class
#' V_data <- NULL # Additional regression covariates
#' glm_form <- "~ 1"
#' 
#' # Run swolca
#' res <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
#'               cluster_id = cluster_id, stratum_id = stratum_id, V_data = V_data,
#'               run_sampler = "both", glm_form = glm_form, adapt_seed = 1,
#'               n_runs = 50, burn = 25, thin = 1, save_res = FALSE)
#'        
#' # Apply variance adjustment to posterior estimates
#' res_adjust <- swolca_var_adjust(res = res, num_reps = 100, save_res = FALSE, 
#'                                 adjust_seed = 1)                        
#' }
swolca_var_adjust <- function(res, alpha = NULL, eta = NULL, mu0 = NULL, 
                              Sig0 = NULL, num_reps = 100, save_res = TRUE,
                              save_path = NULL, adjust_seed = NULL) {
  
  # Begin runtime tracker
  start_time <- Sys.time()
  
  # Set seed
  if (!is.null(adjust_seed)) {
    set.seed(adjust_seed)
  }
  
  #================= Extract dimensions and catch errors =======================
  # Check object class and estimates
  if (!inherits(res, "swolca")) {
    stop("res must be an object of class `swolca`, resulting from a call to the 
         `swolca()` function that includes results from the fixed sampler")
  } else if (is.null(res$estimates)) {
    stop("res must include results from the fixed sampler in the `swolca()` function")
  }
  # Check variance adjustment has not already been performed
  if ("estimates_adjust" %in% names(res)) {
    stop("variance adjustment has already been performed, since res$estimates is not NULL")
  }
  # Check num_reps
  if (!(num_reps %% 1 == 0) | num_reps < 1) {
    stop("num_reps must be a whole number greater than 0, recommended to be at least 50. 
    More replicates will lead to more accurate results but will take longer to run.")
  }
  
  # Extract data elements into the global environment
  K <- res$estimates$K_red
  J <- res$data_vars$J
  R_j <- res$data_vars$R_j
  R <- res$data_vars$R
  n <- res$data_vars$n
  Q <- res$data_vars$Q  # prevent stats::rnorm error
  x_mat <- res$data_vars$x_mat
  y_all <- res$data_vars$y_all
  V <- res$data_vars$V
  w_all <- res$data_vars$w_all
  stratum_id <- res$data_vars$stratum_id
  cluster_id <- res$data_vars$cluster_id
  if (is.null(cluster_id)) {  # no clustering
    cluster_id <- 1:n
  }
  
  # Check hyperparameter dimensions match K
  if (!is.null(alpha)) {
    if (length(alpha) != K) {
      stop("length of alpha must be the same as K")
    }
  }
  if (!is.null(eta)) {
    if ((nrow(eta) != J) | (ncol(eta) != R)) {
      stop("eta must be a matrix with J rows and R columns, where J is
             the number of exposure items and R is the maximum number of 
             exposure categories")
    }
    if (any(eta == 0)) {
      warning("eta has 0 values and may result in rank-difficiency issues
                during the Hessian calculation in the var_adjust() function")
    }
  }
  if (!is.null(mu0)) {
    if (length(mu0) != K | !is.list(mu0) | 
        !(all(lapply(mu0, length) == Q))) {
      stop("mu0 must be a list of length K where each element is a 
           vector of length Q (number of regression covariates excluding latent class)")
    }
  }
  if (!is.null(Sig0)) {
    if (length(Sig0) != K | !is.list(Sig0) | 
        !(all(lapply(Sig0, nrow) == Q)) | 
        !(all(lapply(Sig0, ncol) == Q))) {
      stop("Sig0 must be a list of length K where each element is a 
            QxQ matrix, where Q is the number of regression covariates excluding 
           latent class)")
    }
  }
  
  # Check saving parameters
  if (!is.null(save_res)) {
    if (!is.logical(save_res)) {
      stop("save_res must be a boolean specifying if results should be saved")
    }
    if (save_res) {
      if (is.null(save_path) | !is.character(save_path)) {
        stop("save_path must be a string specifying a path and file name, such as '~/Documents/run'")
      } else {
        last_slash_ind <- regexpr("\\/[^\\/]*$", save_path)
        if (!dir.exists(substr(save_path, start = 1, stop = last_slash_ind))) {
          stop("directory specified in save_path does not exist")
        }
        if (last_slash_ind == length(save_path)) {
          stop("please append the start of a file name to the end of save_path. 
            For example, '~/Documents/run' can produce a saved file named 
            'run_swolca_results.RData'")
        }
      }
    }
  }
  
  #================= Initialize hyperparameters ================================
  # Default hyperparameters for pi and theta
  if (is.null(alpha)) {
    alpha <- rep(1, K) / K   # Hyperparameter for prior for pi
  }
  if (is.null(eta)) {
    # Hyperparameter for prior for theta
    # Unviable categories have value 0.01 to prevent rank deficiency issues
    eta <- matrix(0.01, nrow = J, ncol = R) 
    for (j in 1:J) {
      eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
    }
  }
  # Default hyperparameters for xi
  if (is.null(mu0)) {
    mu0 <- vector("list", K)
    for (k in 1:K) {
      # MVN(0,1) hyperprior for prior mean of xi
      mu0[[k]] <- stats::rnorm(n = Q)
    }
  }
  if (is.null(Sig0)) {
    Sig0 <- vector("list", K)
    for (k in 1:K) {
      # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
      # components and mean variance 2.5 for a weakly informative prior on xi
      Sig0[[k]] <- diag(LaplacesDemon::rinvgamma(n = Q, shape = 3.5, scale = 6.25),
                        nrow = Q, ncol = Q)
    }
  }
  
  #=============== Run Stan model ==============================================
  print("Running variance adjustment")
  
  # Define data for Stan model
  data_stan <- list(K = K, J = J, R = R, n = n, Q = Q, X = x_mat, y = y_all, 
                    V = V, weights = w_all, alpha = alpha, eta = eta, mu0 = mu0, 
                    Sig0 = Sig0)
  
  # Stan parameters of interest
  par_stan <- c('pi', 'theta', 'xi')  # subset of parameters interested in
  
  # Stan model
  mod_stan <- stanmodels$SWOLCA_main

  # Run Stan model
  # Stan will pass warnings from calling 0 chains, but will still create an 
  # out_stan object for the 'grad_log_prob()' method
  out_stan <- rstan::sampling(object = mod_stan, data = data_stan, 
                              pars = par_stan, chains = 0, iter = 0, refresh = 0)
  
  #=============== Convert to unconstrained parameters =========================
  # Convert params from constrained space to unconstrained space
  unc_par_hat <- rstan::unconstrain_pars(out_stan, 
                                         list("pi" = res$estimates$pi_med,
                                              "theta" = res$estimates$theta_med,
                                              "xi" = res$estimates$xi_med))
  # Get posterior MCMC samples in unconstrained space for all parameters
  M <- dim(res$estimates$pi_red)[1]
  unc_par_samps <- lapply(1:M, unconstrain, stan_model = out_stan, K = K, 
                          pi = res$estimates$pi_red, 
                          theta = res$estimates$theta_red, 
                          xi = res$estimates$xi_red)
  unc_par_samps <- matrix(unlist(unc_par_samps), byrow = TRUE, nrow = M)
  
  #=============== Post-processing adjustment in unconstrained space ===========
  # Estimate Hessian
  H_hat <- -1*stats::optimHess(unc_par_hat, 
                               gr = function(x){rstan::grad_log_prob(out_stan, x)})
  
  # Create survey design
  if (!is.null(stratum_id)) {  # Include stratifying variable
    # Survey data frame for specifying survey design
    svy_data <- data.frame(stratum_id = stratum_id, cluster_id = cluster_id,
                           x_mat = x_mat, y_all = y_all, w_all = w_all)
    # Specify survey design
    svydes <- survey::svydesign(ids = ~cluster_id, strata = ~factor(stratum_id), 
                                weights = ~w_all, data = svy_data)
  } else { # No stratifying variable
    # Survey data frame for specifying survey design
    svy_data <- data.frame(cluster_id = cluster_id, x_mat = x_mat, 
                           y_all = y_all, w_all = w_all)
    # Specify survey design
    svydes <- survey::svydesign(ids = ~cluster_id, weights = ~w_all, 
                                data = svy_data)    
  }
  
  # Create svrepdesign
  svyrep <- survey::as.svrepdesign(design = svydes, type = "mrbbootstrap", 
                                   replicates = num_reps)
  # Get survey replicates
  rep_temp <- survey::withReplicates(design = svyrep, theta = grad_par, 
                                     stan_mod = mod_stan, stan_data = data_stan, 
                                     par_stan = par_stan, u_pars = unc_par_hat)
  
  # Compute adjustment
  J_hat <- stats::vcov(rep_temp)
  H_inv <- solve(H_hat)
  V1 <- H_inv %*% J_hat %*% H_inv
  # Check for errors in the matrix inversion
  if (any(is.na(V1))) {
    stop(paste0("NaNs created during variance adjustment, likely due to lack of ",
      "smoothness in the posterior. Please run the sampler for more iterations ", 
      "or do not run the variance adjustment."))
  }
  
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
    V1_pd <- Matrix::nearPD(V1)
    R1 <- chol(V1_pd$mat)
    V1_pd_diff <- sum(abs(eigen(V1)$values - eigen(V1_pd$mat)$values))
    print(paste0("V1: absolute eigenvalue difference to nearest p.d. matrix: ", 
                 V1_pd_diff))
  } else {
    R1 <- chol(V1)
  }
  if (min(Re(eigen(H_inv)$values)) < 0) {
    H_inv_pd <- Matrix::nearPD(H_inv)
    R2_inv <- chol(H_inv_pd$mat)
    H_inv_pd_diff <- sum(abs(eigen(H_inv)$values - eigen(H_inv_pd$mat)$values))
    print(paste0("H_inv: absolute eigenvalue difference to nearest p.d. matrix: ", 
                 H_inv_pd_diff))
    if (H_inv_pd_diff > 5) {
      stop("NaNs created during variance adjustment, likely due to lack of 
      smoothness in the posterior. Please run the sampler for more iterations or 
      do not run the variance adjustment.")
    }
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
  theta_red_adj <- array(NA, dim=c(M, J, K, R))
  xi_red_adj <- array(NA, dim=c(M, K, Q))
  for (i in 1:M) {
    ##### FIX WITH CUSTOMIZED ERROR
    constr_pars <- rstan::constrain_pars(out_stan, par_adj[i,])
    pi_red_adj[i, ] <- constr_pars$pi
    theta_red_adj[i,,,] <- constr_pars$theta
    xi_red_adj[i,,] <- constr_pars$xi
  }
  
  #=============== Output adjusted parameters ==================================
  # Re-normalize pi and theta for each iteration
  pi_red_adj = pi_red_adj / rowSums(pi_red_adj)  
  theta_red_adj <- plyr::aaply(theta_red_adj, c(1, 2, 3), function(x) x / sum(x),
                               .drop = FALSE) 
  
  # Get posterior median estimates
  pi_med_adj <- apply(pi_red_adj, 2, stats::median)
  theta_med_adj <- apply(theta_red_adj, c(2,3,4), stats::median)
  xi_med_adj <- apply(xi_red_adj, c(2,3), stats::median)
  
  # Renormalize posterior median estimates for pi and theta to sum to 1
  pi_med_adj <- pi_med_adj / sum(pi_med_adj)  
  theta_med_adj <- plyr::aaply(theta_med_adj, c(1, 2), function(x) x / sum(x),
                               .drop = FALSE)  # Re-normalize
  
  # Update Phi_med using adjusted xi estimate
  c_all <- res$estimates$c_all
  Phi_med_all_c <- stats::pnorm(V %*% t(xi_med_adj))  # Outcome probabilities for all classes
  Phi_med_adj <- numeric(n)                    # Initialize individual outcome probabilities
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    # Calculate outcome probabilities P(Y=1|-) using updated class assignment
    Phi_med_adj[i] <- Phi_med_all_c[i, c_all[i]] 
  }
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  # Add variance adjustment runtime to overall runtime
  sum_runtime <- runtime + res$runtime
  res$runtime <- sum_runtime
  
  estimates_adjust <- list(pi_red = pi_red_adj, theta_red = theta_red_adj, 
                          xi_red = xi_red_adj, pi_med = pi_med_adj, 
                          theta_med = theta_med_adj, xi_med = xi_med_adj, 
                          Phi_med = Phi_med_adj, c_all = c_all,
                          pred_class_probs = res$estimates$pred_class_probs,
                          log_lik_med = res$estimates$loglik_med)
  
  res$estimates_adjust <- estimates_adjust
  class(res) <- "swolca"
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_swolca_results.RData"))
  }
  
  return(res)
}


