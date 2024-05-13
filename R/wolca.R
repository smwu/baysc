#' Run the WOLCA model
#'
#' @description
#' `wolca` runs an unsupervised weighted overfitted latent class analysis 
#' (WOLCA) described in Wu et al. (2023) and Stephenson et al. (2023), then 
#' saves and returns the results.
#'
#' @inheritParams swolca
#' 
#' @details 
#' If no survey sample adjustments are desired, leave `sampling_wt`, `stratum_id`, 
#' and `cluster_id` to their default `NULL` values.
#' 
#' By default, the function will run two samplers: the adaptive sampler, followed 
#' by the fixed sampler. The adaptive sampler determines the number of latent 
#' classes, which is then used in the fixed sampler for parameter estimation. 
#' If the number of latent classes is already known and only the fixed sampler 
#' needs to be run, specify `"fixed"` for the `run_sampler` argument and specify a 
#' number for `K_fixed`. If only the adaptive sampler is to be run, specify 
#' `"adapt"` for the `run_sampler` argument. Use `adapt_seed` (default is `NULL`) 
#' to specify a seed for the adaptive sampler, and use `fixed_seed` (default is 
#' `NULL`) to specify a separate seed for the fixed sampler. 
#' 
#' `x_mat` is an nxJ matrix with each row corresponding to the J-dimensional 
#' categorical exposure for an individual. If there is no clustering present, 
#' `cluster_id` should be set to the individual IDs. 
#' `K_max` is the maximum number of latent classes allowable, to 
#' be used for the overfitted latent class model if the adaptive sampler is run. 
#' `class_cutoff` is the minimum size of each class as a proportion of the 
#' population, used when determining the number of latent classes.  
#' 
#' To save results, set `save_res = TRUE` (default) and `save_path` to a string
#' that specifies both the location and the beginning of the file name 
#' (e.g., "~/Documents/run"). The file name will have "_wolca_adapt.RData" or 
#' "_wolca_results.RData" appended to it.
#'
#' If hyperparameters for the adaptive or fixed sampler are left as `NULL` 
#' (default), the following default values are used. Let \eqn{K} refer to 
#' `K_max` for the adaptive sampler and `K_fixed` for the fixed sampler. 
#' For \eqn{\pi}, a Dirichlet prior with hyperparameter \eqn{\alpha = 1/K} for 
#' each component. For \eqn{\theta_{jk\cdot}}, a Dirichlet prior with 
#' hyperparameter  \eqn{\eta_j} equal to `rep(1, R_j)` where `R_j` is the number 
#' of categories for exposure item j. If `R_j < R`, the remaining categories have
#' hyperparameter set to 0.01. This is done independently for each exposure item j
#' and is assumed to be the same across latent classes. Note that hyperparameters
#' for the fixed sampler should probably only be specified if running the 
#' fixed sampler directly, bypassing the adaptive sampler. 
#' 
#' To prevent underflow issues, all \eqn{\theta} parameters are restricted to 
#' have a minimum value of \eqn{1e-8}. 
#'
#' @return
#' If the fixed sampler is run, returns an object `res` of class `"wolca"`; a 
#' list containing the following:
#' \describe{
#'   \item{\code{estimates}}{List of posterior model results, resulting from a 
#'   call to [get_estimates_wolca()]}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{List of data variables used, including:
#'   `n`: Sample size.
#'   `J`: Number of exposure items.
#'   `R_j`: Number vector of number of exposure categories for each item; Jx1.
#'   `R`: Maximum number of exposure categories across items.
#'   `w_all`: Vector of sampling weights normalized to sum to n; nx1.
#'   `sampling_wt`: Vector of survey sampling weights; nx1.
#'   `x_mat`: Matrix of multivariate categorical exposures; nxJ.
#'   `stratum_id`: Vector of individual stratum IDs; nx1 or NULL.
#'   `cluster_id`: Vector of individual cluster IDs; nx1 or NULL. 
#'   `ci_level`: Confidence interval level.
#'   }
#'   \item{\code{MCMC_out}}{List of full MCMC output, resulting from a call to 
#'   [run_MCMC_Rcpp()]}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling, resulting 
#'   from a call to [post_process()]}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler}
#' }
#' If `save_res = TRUE` (default), also saves `res` as 
#' `[save_path]_wolca_results.RData`. 
#' 
#' If only the adaptive sampler is run (i.e., `run_sampler` = `"adapt"`), returns
#' list `res` containing:
#' \describe{
#'   \item{\code{MCMC_out}}{List of full MCMC output, resulting from a call to 
#'   [run_MCMC_Rcpp()]}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler, 
#' obtained from the adaptive sampler}
#'   \item{\code{K_MCMC}}{Adaptive sampler MCMC output for K; Mx1, where M is 
#'   the number of MCMC iterations after burn-in and thinning.}
#' }
#' If `save_res = TRUE` (default), also saves `res` as 
#' `[save_path]_wolca_adapt.RData`. 
#' 
#' @importFrom RcppTN rtn
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats rnorm median confint
#' @importFrom survey svydesign svyglm
#' @export
#' 
#' @references 
#' Stephenson, B. J. K., Wu, S. M., Dominici, F. (2023). Identifying dietary 
#' consumption patterns from survey data: a Bayesian nonparametric latent class 
#' model. Journal of the Royal Statistical Society Series A: Statistics in 
#' Society, qnad135.
#' 
#' Wu, S. M., Williams, M. R., Savitsky, T. D., & Stephenson, B. J. 
#' (2023). Derivation of outcome-dependent dietary patterns for low-income women 
#' obtained from survey data using a Supervised Weighted Overfitted Latent Class 
#' Analysis. arXiv preprint arXiv:2310.01575.
#'
#' @examples
#' \dontrun{   
#' # Load data and obtain relevant variables
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
#' }
wolca <- function(x_mat, sampling_wt = NULL, cluster_id = NULL, stratum_id = NULL,  
                  run_sampler = "both", K_max = 30, adapt_seed = NULL, 
                  K_fixed = NULL, fixed_seed = NULL, class_cutoff = 0.05,
                  n_runs = 20000, burn = 10000, thin = 5, update = 10,
                  save_res = TRUE, save_path = NULL,
                  alpha_adapt = NULL, eta_adapt = NULL,
                  alpha_fixed = NULL, eta_fixed = NULL) {
  
  # Begin runtime tracker
  start_time <- Sys.time()

  #================= Read in data ==============================================
  print("Read in data")

  # Obtain dimensions
  n <- dim(x_mat)[1]        # Number of individuals
  J <- dim(x_mat)[2]        # Number of exposure items
  R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
               function(x) length(unique(x)))  
  R <- max(R_j)             # Maximum number of exposure categories across items
  
  # If no sampling weights, set all weights to 1 
  if (is.null(sampling_wt)) {
    w_all <- rep(1, n)
  # Otherwise, obtain normalized weights
  } else {
    kappa <- sum(sampling_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
    w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
  }

  
  #================= Catch errors ==============================================
  catch_errors(x_mat = x_mat, sampling_wt = sampling_wt, cluster_id = cluster_id, 
               stratum_id = stratum_id, run_sampler = run_sampler, K_max = K_max, 
               class_cutoff = class_cutoff, adapt_seed = adapt_seed, 
               fixed_seed = fixed_seed, alpha_adapt = alpha_adapt, 
               eta_adapt = eta_adapt, K_fixed = K_fixed, 
               alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
               n_runs = n_runs, burn = burn, thin = thin, update = update,
               save_res = save_res, save_path = save_path)

  #================= ADAPTIVE SAMPLER ==========================================
  if (run_sampler %in% c("both", "adapt")) { # Run adaptive sampler
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
      # Hyperparameter for prior for theta
      # Unviable categories have value 0.01 to prevent rank deficiency issues
      eta_adapt <- matrix(0.01, nrow = J, ncol = R) 
      for (j in 1:J) {
        eta_adapt[j, 1:R_j[j]] <- rep(1, R_j[j]) 
      }
    }
    
    #================= Initialize OLCA model =====================================
    # Obtain pi, theta, c_all
    OLCA_params <- init_OLCA(alpha = alpha_adapt, eta = eta_adapt, n = n,
                             K = K_max, J = J, R = R)
    
    #================= Run adaptive sampler to obtain number of classes ========
    # Obtain pi_MCMC, theta_MCMC, c_all_MCMC
    MCMC_out <- run_MCMC_Rcpp_wolca(OLCA_params = OLCA_params, n_runs = n_runs, 
                                    burn = burn, thin = thin, K = K_max, J = J, 
                                    R = R, n = n, w_all = w_all, x_mat = x_mat, 
                                    alpha = alpha_adapt, eta = eta_adapt, 
                                    update = update)
    #================= Post-processing for adaptive sampler ====================
    # Get median number of classes with >= cutoff% of individuals, over all iterations
    M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
    K_MCMC <- rowSums(MCMC_out$pi_MCMC >= class_cutoff)
    K_med <- round(stats::median(K_MCMC))
    # Get number of unique classes for fixed sampler
    K_fixed <- K_med
    print(paste0("K_fixed: ", K_fixed))
    
    # Create adaptive output list (fixed sampler replaces this if run)
    res <- list(MCMC_out = MCMC_out, K_fixed = K_fixed, K_MCMC = K_MCMC)
    # Save output
    if (save_res) {
      save(res, file = paste0(save_path, "_wolca_adapt.RData"))
    }
    
    # Reduce memory burden
    rm(OLCA_params, MCMC_out)
  }  

  #================= FIXED SAMPLER =============================================
  if (run_sampler %in% c("both", "fixed")) {
    print("Running fixed sampler...")
    
    # Catch errors: check hyperparameter dimensions for fixed sampler
    catch_errors(x_mat = x_mat, K_fixed = K_fixed, 
                 alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
                 n_runs = n_runs, burn = burn, thin = thin, update = update)
    
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
      # Hyperparameter for prior for theta
      # Unviable categories have value 0.01 to prevent rank deficiency issues
      eta_fixed <- matrix(0.01, nrow = J, ncol = R) 
      for (j in 1:J) {
        eta_fixed[j, 1:R_j[j]] <- rep(1, R_j[j]) 
      }
    }
    
    #================= Run fixed sampler to obtain posteriors ====================
    # Initialize OLCA model using fixed number of classes. Obtain pi, theta, c_all
    OLCA_params <- init_OLCA(alpha = alpha_fixed, eta = eta_fixed, n = n,
                             K = K_fixed, J = J, R = R)
  
    # Run MCMC algorithm using fixed number of classes
    # Obtain pi_MCMC, theta_MCMC, c_all_MCMC
    MCMC_out <- run_MCMC_Rcpp_wolca(OLCA_params = OLCA_params, n_runs = n_runs, 
                                    burn = burn, thin = thin, K = K_fixed, J = J, 
                                    R = R, n = n, w_all = w_all, x_mat = x_mat, 
                                    alpha = alpha_fixed, eta = eta_fixed,
                                    update = update)
  
    # Post-processing to recalibrate labels and remove extraneous empty classes
    # Obtain K_med, pi, theta
    post_MCMC_out <- post_process_wolca(MCMC_out = MCMC_out, J = J, R = R,
                                        class_cutoff = class_cutoff)
  
    # Obtain posterior estimates, reduce number of classes, analyze results
    # Obtain K_red, pi_red, theta_red, pi_med, theta_med, c_all, pred_class_probs
    estimates <- get_estimates_wolca(MCMC_out = MCMC_out, 
                                     post_MCMC_out = post_MCMC_out, n = n, J = J,
                                     x_mat = x_mat)
    
    # Create output list. Replaces adaptive sampler output list
    res <- list(estimates = estimates, MCMC_out = MCMC_out,
                post_MCMC_out = post_MCMC_out, K_fixed = K_fixed)
  }
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  res$runtime <- runtime
  
  # Store data variables used
  data_vars <- list(n = n, J = J, R_j = R_j, R = R, w_all = w_all, 
                    sampling_wt = sampling_wt, x_mat = x_mat,
                    stratum_id = stratum_id, cluster_id = cluster_id)
  res$data_vars <- data_vars
  
  class(res) <- "wolca"
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_wolca_results.RData"))
  }
  
  return(res)
}




