#' Run the WOLCA model
#'
#' @description
#' `wolca` runs a two-step model with an unsupervised weighted overfitted latent 
#' class analysis (WOLCA) in the first step and saves and returns the results.
#'
#' @inheritParams swolca
#' @param glm_form String specifying formula for probit regression, 
#' including latent class as a covariate. For example, `"~ c_all"` for the model 
#' with only latent class as covariates. All variables in `glm_form` must be 
#' found in `V`.
#' 
#' @details 
#' `wolca` is a two-step approach that runs an unsupervised WOLCA in the first
#' step to derive latent class patterns and subsequently treats the class 
#' assignments as fixed and includes them as covariates in a frequentist 
#' survey-weighted probit regression model that uses an asymptotic sandwich 
#' estimator for variance estimation. 
#' 
#' By default, the function will run both samplers for the WOLCA step, running 
#' the adaptive sampler first to determine the number of latent classes, and 
#' then using the determined number of latent classes to run the fixed sampler 
#' for parameter estimation of the unsupervised WOLCA, followed by the 
#' survey-weighted regression model. If the number of latent classes is already 
#' known and only the fixed sampler is to be run, specify `"fixed"` for the 
#' `run_sampler` argument and specify a number for `K_fixed`. Id only the 
#' adaptive sampler is to be run, specify `"adapt"` for the `run_sampler` 
#' argument. Use `adapt_seed` (default is `NULL`) to specify a seed 
#' for the adaptive sampler, and use `fixed_seed` (default is `NULL`) to specify 
#' a separate seed for the fixed sampler.
#' 
#' `x_mat` is an nxp matrix with each row corresponding to the J-dimensional 
#' categorical exposure for an individual. If there is no clustering present, 
#' `cluster_id` should be set to the individual IDs. `V` is the design matrix for 
#' the probit regression, including the intercept and all covariates other than 
#' latent class. `K_max` is the maximum number of latent classes allowable, to 
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
#' hyperparameter \eqn{\eta = 1} for each component. Note that hyperparameters
#' for the fixed sampler should probably only be specified if running the 
#' fixed sampler directly, bypassing the adaptive sampler. 
#'
#' @return
#' If the fixed sampler is run, returns list `res` containing:
#' \describe{
#'   \item{\code{estimates}}{List of posterior model results}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{List of data variables used}
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler}
#'   \item{\code{K_MCMC}}{If `K_fixed = NULL` and the adaptive sampler is run,
#'   output list also contains MCMC output for the number of classes with size
#'   greater than `class_cutoff` for each iteration}
#' }
#' If `save_res = TRUE` (default), also saves `res` as 
#' `[save_path]_wolca_results.RData`. 
#' 
#' If only the adaptive sampler is run (i.e., `run_sampler` = `"adapt"`), returns
#' list `res` containing:
#' \describe{
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler, 
#' obtained from the adaptive sampler}
#'   \item{\code{K_MCMC}}{Adaptive sampler MCMC output for K}
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
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
#' stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
#' sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
#' n <- dim(x_mat)[1]                   # Number of individuals
#' 
#' # Probit model only includes latent class
#' V <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
#' # Survey-weighted regression formula
#' glm_form <- "~ c_all"
#' 
#' # Run wolca
#' res <- wolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
#'        cluster_id = cluster_id, stratum_id = stratum_id, V = V, 
#'        run_sampler = "both", glm_form = glm_form, adapt_seed = 1, 
#'        n_runs = 50, burn = 25, thin = 1, save_res = FALSE)
#'
wolca <- function(x_mat, y_all, sampling_wt, cluster_id, stratum_id, 
                  V, run_sampler = "both", glm_form, K_max = 30, 
                  adapt_seed = NULL, class_cutoff = 0.05,
                  alpha_adapt = NULL, eta_adapt = NULL,
                  alpha_fixed = NULL, eta_fixed = NULL,
                  K_fixed = NULL, fixed_seed = NULL,
                  n_runs = 20000, burn = 10000, thin = 5, 
                  save_res = TRUE, save_path = NULL) {
  
  # Begin runtime tracker
  start_time <- Sys.time()

  #================= Read in data ==============================================
  print("Read in data")

  # Obtain dimensions
  n <- dim(x_mat)[1]        # Number of individuals
  J <- dim(x_mat)[2]        # Number of exposure items
  R <- max(apply(x_mat, 2,  # Number of exposure categories
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM

  # Obtain normalized weights
  kappa <- sum(sampling_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
  w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
  
  #================= Catch errors ==============================================
  catch_errors(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
               cluster_id = cluster_id, stratum_id = stratum_id, V = V,
               run_sampler = run_sampler, glm_form = glm_form,
               K_max = K_max, class_cutoff = class_cutoff,
               alpha_adapt = alpha_adapt, eta_adapt = eta_adapt, 
               K_fixed = K_fixed, alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
               n_runs = n_runs, burn = burn, thin = thin, 
               save_res = save_res, save_path = save_path)
  
  q <- ncol(V)  # Number of regression covariates excluding class assignment   
  
  if (!grepl("c_all", glm_form)) {
    stop("glm_form should include latent class assignments, c_all")
  }

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
      eta_adapt <- rep(1, R)                 # Hyperparameter for prior for theta
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
                                    alpha = alpha_adapt, eta = eta_adapt)
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
      save(res, file = paste0(save_path, "_swolca_adapt.RData"))
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
                 n_runs = n_runs, burn = burn, thin = thin)
    
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
      eta_fixed <- rep(1, R)                 # Hyperparameter for prior for theta
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
                                    alpha = alpha_fixed, eta = eta_fixed)
  
    # Post-processing to recalibrate labels and remove extraneous empty classes
    # Obtain K_med, pi, theta
    post_MCMC_out <- post_process_wolca(MCMC_out = MCMC_out, J = J, R = R,
                                        class_cutoff = class_cutoff)
  
    # Obtain posterior estimates, reduce number of classes, analyze results
    # Obtain K_red, pi_red, theta_red, pi_med, theta_med, c_all, pred_class_probs
    estimates <- get_estimates_wolca(MCMC_out = MCMC_out, 
                                     post_MCMC_out = post_MCMC_out, n = n, J = J,
                                     x_mat = x_mat)
  
    #================= Fit probit model ==========================================
    
    estimates <- fit_probit_wolca(estimates = estimates, glm_form = glm_form, 
                                  stratum_id = stratum_id, cluster_id = cluster_id, 
                                  x_mat = x_mat, y_all = y_all, w_all = w_all, 
                                  V = V, q = q)
    
    # Create output list. Replaces adaptive sampler output list
    res <- list(estimates = estimates, V = V, MCMC_out = MCMC_out,
                post_MCMC_out = post_MCMC_out, K_fixed = K_fixed)
  }
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  res$runtime <- runtime
  
  # Store data variables used
  data_vars <- list(n = n, J = J, R = R, q = q, sample_wt = sampling_wt,
                    X_data = x_mat, Y_data = y_all, V = V,
                    true_Si = stratum_id, cluster_id = cluster_id)
  res$data_vars <- data_vars
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_wolca_results.RData"))
  }
  
  return(res)
}




