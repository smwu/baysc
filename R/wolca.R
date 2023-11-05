#' Run the WOLCA model
#'
#' @description
#' `wolca` runs a two-step model with an unsupervised weighted overfitted latent 
#' class analysis (WOLCA) in the first step and saves and returns the results.
#'
#' @inheritParams swolca
#' @param glm_form String specifying formula to use for probit regression. For 
#' example, `"y_all ~ c_all"` for the model with only latent class as a covariate. 
#' 
#' @details 
#' If hyperparameters are left as `NULL` (default), the following values are 
#' used. 
#'
#' @return
#' Returns list `res` containing:
#' \describe{
#'   \item{\code{estimates}}{List of posterior model results}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{List of data variables used}
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler}
#'   \item{\code{K_MCMC}}{If `K_true = NULL` and the adaptive sampler is run,
#'   output list also contains MCMC output for the number of classes with size
#'   greater than `class_cutoff` for each iteration}
#' }
#'
#' If `save_res = TRUE` (default), also saves `res` as `[save_path]_wolca_results.RData`
#' and, if `K_true = NULL` so that the adaptive sampler is run, list `adapt_MCMC`
#' is saved as  `[save_path]_wolca_adapt.RData`. List `adapt_MCMC` contains:
#' \describe{
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler, 
#' obtained from the adaptive sampler}
#'   \item{\code{K_MCMC}}{Adaptive sampler MCMC output for K}
#' }
#' 
#' @importFrom RcppTN rtn
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats rnorm median confint
#' @importFrom survey svydesign svyglm
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
#' stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
#' sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
#' n <- dim(x_mat)[1]        # Number of individuals
#' 
#' # Probit model only includes latent class
#' V <- matrix(1, nrow = n) # Regression design matrix without class assignment
#' # Survey-weighted regression formula
#' glm_form <- "y_all ~ c_all"
#'
#' # Stan model
#' mod_stan <- stanmodels$SWOLCA_main
#' 
#' # Run swolca
#' res <- wolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
#'        cluster_id = cluster_id, stratum_id = stratum_id, V = V, 
#'        glm_form = glm_form, adapt_seed = 1, n_runs = 50, burn = 25, thin = 1, 
#'        mod_stan = mod_stan, save_res = FALSE)
#'
wolca <- function(x_mat, y_all, sampling_wt, cluster_id, stratum_id, 
                  V, glm_form, K_max = 30, 
                  adapt_seed = NULL, class_cutoff = 0.05,
                  alpha_adapt = NULL, eta_adapt = NULL,
                  mu0_adapt = NULL, Sig0_adapt = NULL,
                  alpha_fixed = NULL, eta_fixed = NULL,
                  mu0_fixed = NULL, Sig0_fixed = NULL,
                  K_true = NULL, fixed_seed = NULL,
                  n_runs = 20000, burn = 10000, thin = 5, mod_stan,
                  save_res = TRUE, save_path = NULL) {
  
  # Begin runtime tracker
  start_time <- Sys.time()

  #================= Read in data ==============================================
  print("Read in data")

  # Obtain dimensions
  n <- dim(data_vars$X_data)[1]        # Number of individuals
  p <- dim(data_vars$X_data)[2]        # Number of exposure items
  d <- max(apply(data_vars$X_data, 2,  # Number of exposure categories
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
  q <- ncol(V)              # Number of regression covariates excluding class assignment

  # Obtain normalized weights
  kappa <- sum(data_vars$sample_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
  w_all <- c(data_vars$sample_wt / kappa) # Weights normalized to sum to n, nx1

  #================= ADAPTIVE SAMPLER ==========================================

  if (!is.null(K_true)) {  # True number of classes is known
    
    K_fixed <- K_true
  } else {  # Run adaptive sampler
    print("Running adaptive sampler...")
    
    # Set seed
    if(!is.null(adapt_seed)) {
      set.seed(adapt_seed)
    }
    
    #================= Initialize hyperparameters ==============================
    # Default hyperparameters for pi and theta
    if (is.null(alpha_adapt)) {
      alpha_adapt <- rep(1, K_max) / K_max   # Hyperparameter for prior for pi
    }
    if (is.null(eta_adapt)) {
      eta_adapt <- rep(1, d)                 # Hyperparameter for prior for theta
    }
    
    #================= Initialize OLCA model =====================================
    # Obtain pi, theta, c_all
    OLCA_params <- init_OLCA(alpha = alpha_adapt, eta = eta_adapt, n = n,
                             K = K_max, p = p, d = d)
    
    #================= Run adaptive sampler to obtain number of classes ========
    # Obtain pi_MCMC, theta_MCMC, c_all_MCMC
    MCMC_out <- run_MCMC_Rcpp_wolca(OLCA_params = OLCA_params, n_runs = n_runs, 
                                    burn = burn, thin = thin, K = K_max, p = p, 
                                    d = d, n = n, w_all = w_all, x_mat = x_mat, 
                                    alpha = alpha_adapt, eta = eta_adapt)
    #================= Post-processing for adaptive sampler ====================
    # Get median number of classes with >= cutoff% of individuals, over all iterations
    M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
    K_MCMC <- rowSums(MCMC_out$pi_MCMC >= class_cutoff)
    K_med <- round(stats::median(K_MCMC))
    # Get number of unique classes for fixed sampler
    K_fixed <- K_med
    print(paste0("K_fixed: ", K_fixed))
    # Save adaptive output
    adapt_MCMC <- list(MCMC_out = MCMC_out, K_fixed = K_fixed, K_MCMC = K_MCMC)
    if (save_res) {
      save(adapt_MCMC, file = paste0(save_path, "_wolca_adapt.RData"))
    }
    # Reduce memory burden
    rm(OLCA_params, MCMC_out)
  }  

  #================= FIXED SAMPLER =============================================
  print("Running fixed sampler...")
  
  # Set seed
  if (!is.null(fixed_seed)) {
    set.seed(fixed_seed)
  }
  
  #================= Initialize hyperparameters ==============================
  # Default hyperparameters for pi and theta
  if (is.null(alpha_fixed)) {
    alpha_fixed <- rep(1, K_fixed) / K_fixed   # Hyperparameter for prior for pi
  }
  if (is.null(eta_fixed)) {
    eta_fixed <- rep(1, d)                 # Hyperparameter for prior for theta
  }
  
  #================= Run fixed sampler to obtain posteriors ====================
  # Initialize OLCA model using fixed number of classes. Obtain pi, theta, c_all
  OLCA_params <- init_OLCA(alpha = alpha_fixed, eta = eta_fixed, n = n,
                           K = K_fixed, p = p, d = d)

  # Run MCMC algorithm using fixed number of classes
  # Obtain pi_MCMC, theta_MCMC, c_all_MCMC
  MCMC_out <- run_MCMC_Rcpp_wolca(OLCA_params = OLCA_params, n_runs = n_runs, 
                                  burn = burn, thin = thin, K = K_fixed, p = p, 
                                  d = d, n = n, w_all = w_all, x_mat = x_mat, 
                                  alpha = alpha_fixed, eta = eta_fixed)

  # Post-processing to recalibrate labels and remove extraneous empty classes
  # Obtain K_med, pi, theta
  post_MCMC_out <- post_process_wolca(MCMC_out = MCMC_out, p = p, d = d)

  # Obtain posterior estimates, reduce number of classes, analyze results
  # Obtain K_red, pi_red, theta_red, pi_med, theta_med, c_all, pred_class_probs
  estimates <- get_estimates_wolca(MCMC_out = MCMC_out, 
                                   post_MCMC_out = post_MCMC_out, n = n, p = p,
                                   x_mat = x_mat)

  #================= Fit probit model ==========================================
  
  estimates <- fit_probit_wolca(estimates = estimates, glm_form = glm_form, 
                                stratum_id = stratum_id, cluster_id = cluster_id, 
                                x_mat = x_mat, y_all = y_all, w_all = w_all, 
                                V = V, q = q)

  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  
  # Store data variables used
  data_vars <- list(n = n, p = p, d = d, q = q, sample_wt = sampling_wt,
                    X_data = x_mat, Y_data = y_all, V = V,
                    true_Si = stratum_id, cluster_id = cluster_id)
  # Create output list
  res <- list(estimates = estimates, runtime = runtime,
              data_vars = data_vars, V = V, MCMC_out = MCMC_out,
              post_MCMC_out = post_MCMC_out, K_fixed = K_fixed)
  if (is.null(K_true)) {
    res$K_MCMC <- K_MCMC
  }
  if (save_res) {
    save(res, file = paste0(save_path, "_wolca_results.RData"))
  }
  
  return(res)
}




