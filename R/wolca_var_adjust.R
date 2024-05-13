#' Unconstrain parameters for WOLCA
#'
#' `unconstrain_wolca` is a helper function that converts each row of an input array 
#' of MCMC parameter output from constrained space to unconstrained space in Stan
#' 
#' @param i Row index
#' @param K Number of classes
#' @param stan_model Stan model
#' @param pi MCMC matrix output for pi; MxK
#' @param theta MCMC array output for theta; MxJxKxR
#' 
#' @return Outputs vector of unconstrained parameters
#' @importFrom rstan unconstrain_pars
#' @keywords internal
#' @export
unconstrain_wolca <- function(i, K, stan_model, pi, theta) {
  u_pars <- rstan::unconstrain_pars(stan_model, 
                                    list("pi" = pi[i,], "theta" = theta[i,,,]))
  return(u_pars)
}



#' Post-processing variance adjustment
#' 
#' @description 
#' `wolca_var_adj` applies applies the post-processing variance adjustment after 
#' a call to [wolca()] to correct for underestimation of posterior intervals.
#' 
#' @inheritParams swolca_var_adjust
#' @param res An object of class `"wolca"`, resulting from a call to [wolca()],
#' containing the unadjusted estimates.
#' @param save_path String specifying directory and file name to save results, 
#' e.g., "~/Documents/run". Default is `NULL`. If this is the same as the file 
#' name specified in [wolca()], the unadjusted results are overwritten with the 
#' adjusted results.
#' 
#' @details
#' `wolca_var_adjust` applies a post-processing variance adjustment that rescales 
#' the variance to obtain correct coverage of posterior intervals, adapted from 
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
#' (e.g., "~/Documents/run"). The file name will have "_wolca_results.RData" 
#' appended to it, overwriting the unadjusted results if the file names are the 
#' same.
#' 
#' If hyperparameters are left as `NULL` (default), the following default 
#' values are used. Let \eqn{K} refer to the final number of latent class 
#' obtained from running [wolca()], available at `res$estimates$K_red`.
#' For \eqn{\pi}, a Dirichlet prior with hyperparameter \eqn{\alpha = 1/K} for 
#' each component. For \eqn{\theta_{jk\cdot}}, a Dirichlet prior with 
#' hyperparameter  \eqn{\eta_j} equal to `rep(1, R_j)` where `R_j` is the number 
#' of levels for exposure item j. If `R_j < R`, the remaining levels have
#' hyperparameter set to 0.01. This is done independently for each exposure item j
#' and is assumed to be the same across latent classes. 
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
#' Returns an object `res` of class `"wolca"`, which includes all outputs from 
#' [wolca()] as well as a list `estimates_adjust` containing:
#' \describe{
#'   \item{\code{pi_red}}{Matrix of adjusted posterior samples for pi. Mx(K_red), 
#'   where M is the number of MCMC iterations after burn-in and thinning.}
#'   \item{\code{theta_red}}{Array of adjusted posterior samples for theta. MxJx(K_red)xR}
#'   \item{\code{pi_med}}{Vector of adjusted posterior median estimates for pi. (K_red)x1}
#'   \item{\code{theta_med}}{Array of adjusted posterior median estimates for theta. px(K_red)xR}
#'   \item{\code{c_all}}{Vector of final individual class assignments from `wolca()`. nx1}
#'   \item{\code{pred_class_probs}}{Matrix of individual posterior class 
#'   probabilities from `wolca()`. nx(K_red)}
#' }
#' The `runtime` output for `res` is also updated to include the runtime for the 
#' variance adjustment in addition to the runtime for the main `wolca()` model.
#' 
#' If `save_res = TRUE` (default), the updated `res` object is saved as 
#' `[save_path]_wolca_results.RData`, overwriting the unadjusted results if the 
#' file names are the same. 
#'
#' @seealso [wolca()]
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
#' cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
#' stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
#' sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
#' n <- dim(x_mat)[1]                   # Number of individuals
#' 
#' # Run wolca
#' res <- wolca(x_mat = x_mat, sampling_wt = sampling_wt, 
#'              cluster_id = cluster_id, stratum_id = stratum_id, 
#'              run_sampler = "both", adapt_seed = 1, n_runs = 50, burn = 25, 
#'              thin = 1, save_res = FALSE)
#'        
#' # Apply variance adjustment to posterior estimates
#' res_adjust <- wolca_var_adjust(res = res, num_reps = 100, save_res = FALSE, 
#'                                adjust_seed = 1)                        
#' }
wolca_var_adjust <- function(res, alpha = NULL, eta = NULL, num_reps = 100, 
                             save_res = TRUE, save_path = NULL, 
                             adjust_seed = NULL) {
  
  # Begin runtime tracker
  start_time <- Sys.time()
  
  # Set seed
  if (!is.null(adjust_seed)) {
    set.seed(adjust_seed)
  }
  
  #================= Extract dimensions and catch errors =======================
  # Check object class and estimates
  if (!inherits(res, "wolca")) {
    stop("res must be an object of class `wolca`, resulting from a call to the 
         `wolca()` function that includes results from the fixed sampler")
  } else if (is.null(res$estimates)) {
    stop("res must include results from the fixed sampler in the `wolca()` function")
  }
  # Check variance adjustment has not already been performed
  if ("estimates_adjust" %in% names(res)) {
    stop("variance adjustment has already been performed, since res$estimates is not NULL")
  }
  # Check num_reps
  if ((num_reps %% 100 != 0) | num_reps < 1) {
    stop("num_reps must be a whole number greater than 0, recommended to be at least 50. 
    More replicates will lead to more accurate results but will take longer to run.")
  }
  
  # Extract data elements into the global environment
  K <- res$estimates$K_red
  J <- res$data_vars$J
  R_j <- res$data_vars$R_j
  R <- res$data_vars$R
  n <- res$data_vars$n
  x_mat <- res$data_vars$x_mat
  w_all <- res$data_vars$w_all
  stratum_id <- res$data_vars$stratum_id
  cluster_id <- res$data_vars$cluster_id
  if (is.null(cluster_id)) {  # no clustering
    cluster_id <- 1:n
  }
  
  # Get final number of classes (usually fewer classes than K_fixed)
  K <- res$estimates$K_red
  
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
            'run_wolca_results.RData'")
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
  
  #=============== Run Stan model ==============================================
  print("Running variance adjustment")
  
  # Define data for Stan model
  data_stan <- list(K = K, J = J, R = R, n = n, X = x_mat, weights = w_all, 
                    alpha = alpha, eta = eta)
  
  # Stan parameters of interest
  par_stan <- c('pi', 'theta')  # subset of parameters interested in
  
  # Stan model
  mod_stan <- stanmodels$WOLCA_main
  
  # Run Stan model
  # Stan will pass warnings from calling 0 chains, but will still create an 
  # out_stan object for the 'grad_log_prob()' method
  out_stan <- rstan::sampling(object = mod_stan, data = data_stan, 
                              pars = par_stan, chains = 0, iter = 0, refresh = 0)
  
  #=============== Convert to unconstrained parameters =========================
  # Convert params from constrained space to unconstrained space
  unc_par_hat <- rstan::unconstrain_pars(out_stan, 
                                         list("pi" = res$estimates$pi_med,
                                              "theta" = res$estimates$theta_med))
  # Get posterior MCMC samples in unconstrained space for all parameters
  M <- dim(res$estimates$pi_red)[1]
  unc_par_samps <- lapply(1:M, unconstrain_wolca, stan_model = out_stan, K = K, 
                          pi = res$estimates$pi_red, 
                          theta = res$estimates$theta_red)
  unc_par_samps <- matrix(unlist(unc_par_samps), byrow = TRUE, nrow = M)
  
  #=============== Post-processing adjustment in unconstrained space ===========
  # Estimate Hessian
  H_hat <- -1*stats::optimHess(unc_par_hat, 
                               gr = function(x){rstan::grad_log_prob(out_stan, x)})
  
  # Create survey design
  if (!is.null(stratum_id)) {  # Include stratifying variable
    # Survey data frame for specifying survey design
    svy_data <- data.frame(stratum_id = stratum_id, cluster_id = cluster_id,
                           x_mat = x_mat, w_all = w_all)
    # Specify survey design
    svydes <- survey::svydesign(ids = ~cluster_id, strata = ~factor(stratum_id), 
                                weights = ~w_all, data = svy_data)
  } else { # No stratifying variable
    # Survey data frame for specifying survey design
    svy_data <- data.frame(cluster_id = cluster_id, x_mat = x_mat, w_all = w_all)
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
      do not run variance adjustment.")
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
  for (i in 1:M) {
    ##### FIX WITH CUSTOMIZED ERROR
    constr_pars <- rstan::constrain_pars(out_stan, par_adj[i,])
    pi_red_adj[i, ] <- constr_pars$pi
    theta_red_adj[i,,,] <- constr_pars$theta
  }
  
  #=============== Output adjusted parameters ==================================
  # Re-normalize pi and theta for each iteration
  pi_red_adj = pi_red_adj / rowSums(pi_red_adj)  
  theta_red_adj <- plyr::aaply(theta_red_adj, c(1, 2, 3), function(x) x / sum(x),
                               .drop = FALSE) 
  
  # Get posterior median estimates
  pi_med_adj <- apply(pi_red_adj, 2, stats::median)
  theta_med_adj <- apply(theta_red_adj, c(2,3,4), stats::median)
  
  # Renormalize posterior median estimates for pi and theta to sum to 1
  pi_med_adj <- pi_med_adj / sum(pi_med_adj)  
  theta_med_adj <- plyr::aaply(theta_med_adj, c(1, 2), function(x) x / sum(x),
                               .drop = FALSE)  # Re-normalize
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  # Add variance adjustment runtime to overall runtime
  sum_runtime <- runtime + res$runtime
  res$runtime <- sum_runtime
  
  estimates_adjust <- list(pi_red = pi_red_adj, theta_red = theta_red_adj, 
                          pi_med = pi_med_adj, theta_med = theta_med_adj, 
                          c_all = res$estimates$c_all,
                          pred_class_probs = res$estimates$pred_class_probs)
  
  res$estimates_adjust <- estimates_adjust
  class(res) <- "wolca"
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_wolca_results.RData"))
  }
  
  return(res)
}


