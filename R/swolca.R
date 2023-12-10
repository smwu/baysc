#' Run the SWOLCA model
#'
#' @description
#' `swolca` runs a supervised weighted overfitted latent class analysis (SWOLCA)
#' and saves and returns the results.
#'
#' @param x_mat Matrix of multivariate categorical exposures. nxJ
#' @param y_all Vector of outcomes. nx1
#' @param sampling_wt Vector of survey sampling weights. nx1. If no sampling 
#' weights are available, set this to a vector of ones. 
#' @param cluster_id Vector of individual cluster IDs. nx1. Default is `NULL`,
#' indicating each individual is their own cluster.
#' @param stratum_id Vector of individual stratum IDs. nx1. Default is `NULL`,
#' indicating no stratification.
#' @param V_data Dataframe of additional regression covariates. nxq. Factor 
#' covariates must be converted to factors. If `NULL` (default), no additional 
#' covariates are to be included. All variables in `glm_form` must 
#' be found in `V_data`.
#' @param run_sampler String specifying which sampler(s) should be run. Must be 
#' one of `"both"` (default), `"fixed"`, or `"adapt"`.
#' @param glm_form String specifying formula for probit regression, excluding 
#' outcome and latent class. For example, `"~ 1"` for the model with only 
#' latent class as covariates. All variables in `glm_form` must be found in `V_data`.
#' Do not specify interaction terms for latent class by additional covariates, 
#' as these terms are already included. 
#' @param K_max Upper limit for number of classes. Default is 30.
#' @param adapt_seed Numeric seed for adaptive sampler. Default is `NULL`.
#' @param class_cutoff Minimum class size proportion when determining number of
#' classes in adaptive sampler. Default is 0.05.
#' @param alpha_adapt Adaptive sampler hyperparameter for prior for class 
#' membership probabilities \eqn{\pi}. Default is `NULL` and default values are 
#' used. If specified, must be (`K_max`)x1. 
#' @param eta_adapt Adaptive sampler hyperparameter for prior for item consumption
#' level probabilities \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class 
#' \eqn{k}, assumed to be the same across classes. Default is `NULL` and default 
#' values are used. If specified, must be JxR, where J is the number of exposure
#' items and R is the maximum number of categories for the exposure.
#' @param mu0_adapt Adaptive sampler mean hyperparameters for regression coefficients
#' \eqn{\xi_{k\cdot}} for each class \eqn{k}. Default is `NULL` and default values are used. 
#' If specified, must be a list of `K_max` vectors of dimension qx1, where q is 
#' the number of regression covariates excluding latent class assignment.  
#' @param Sig0_adapt Adaptive sampler variance hyperparameters for regression coefficients
#' \eqn{\xi_{k\cdot}} for each class \eqn{k}. Default is `NULL` and default values are 
#' used. If specified, must be a list of `K_max` qxq matrices, where q is the 
#' number of regression covariates excluding latent class assignment. 
#' @param fixed_seed Numeric seed for fixed sampler. Default is `NULL`.
#' @param K_fixed True number of classes, if known. Default is `NULL`, as it is
#' not necessary for the adaptive sampler. If bypassing the adaptive sampler and
#' running the fixed sampler directly, need to specify a value here.
#' @param alpha_fixed Fixed sampler hyperparameter for prior for \eqn{\pi}. Default is 
#' `NULL` and default values are used. If specified, must be (`K_fixed`)x1. 
#' @param eta_fixed Fixed sampler hyperparameter for prior for item consumption
#' level probabilities \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class 
#' \eqn{k}, assumed to be the same across classes. Default is `NULL` and default 
#' values are used. If specified, must be JxR, where J is the number of exposure
#' items and R is the maximum number of categories for the exposure.
#' @param mu0_fixed Fixed sampler mean hyperparameters for regression coefficients
#' \eqn{\xi_{k\cdot}} for each class \eqn{k}. Default is `NULL` and default values are used. 
#' If specified, must be a list of `K_fixed` vectors of dimension qx1, where q is 
#' the number of regression covariates excluding latent class assignment. 
#' @param Sig0_fixed Fixed sampler variance hyperparameters for regression coefficients
#' \eqn{\xi_{k\cdot}} for each class \eqn{k}. Default is `NULL` and default values are 
#' used. If specified, must be a list of `K_fixed` qxq matrices, where q is the 
#' number of regression covariates excluding latent class assignment. 
#' @param n_runs Number of MCMC iterations. Default is 20000.
#' @param burn Number of MCMC iterations to drop as a burn-in period. Default is 10000.
#' @param thin Thinning factor for MCMC iterations. Default is 5.
#' @param adjust_var Boolean for the fixed sampler specifying if the 
#' post-processing variance adjustment for accurate interval coverage should be 
#' applied. Default = `TRUE`.
#' @param num_reps Number of bootstrap replicates to use for the variance 
#' adjustment estimate. Default is 100.
#' @param save_res Boolean specifying if results should be saved. Default = `TRUE`.
#' @param save_path String specifying directory and file name to save results, 
#' e.g., "~/Documents/run". Default is `NULL`.
#' 
#' @details 
#' 
#' By default, the function will run both samplers, running the adaptive sampler 
#' first to determine the number of latent classes, and then using the determined 
#' number of latent classes to run the fixed sampler for parameter estimation. 
#' If the number of latent classes is already known and only the fixed sampler
#' is to be run, specify `"fixed"` for the `run_sampler` argument and specify a 
#' number for `K_fixed`. Id only the adaptive sampler is to be run, specify 
#' `"adapt"` for the `run_sampler` argument. 
#' Use `adapt_seed` (default is `NULL`) to specify a seed 
#' for the adaptive sampler, and use `fixed_seed` (default is `NULL`) to specify 
#' a separate seed for the fixed sampler. Specify `adjust_var = TRUE` (default)
#' to adjust the estimates to prevent overly conservative interval estimates. 
#' 
#' `x_mat` is an nxJ matrix with each row corresponding to the J-dimensional 
#' categorical exposure for an individual. If there is no clustering present, 
#' `cluster_id` should be set to the individual IDs. `V_data` includes all 
#' additional covariates other than latent class that are to be included in the
#' probit regression model. If there are no additional covariates, `V_data` 
#' should be `NULL` (default). 
#' `K_max` is the maximum number of latent classes allowable, to be used for 
#' the overfitted latent class model if the adaptive sampler is run. 
#' `class_cutoff` is the minimum size of each class as a proportion of the 
#' population, used when determining the number of latent classes.  
#' 
#' To save results, set `save_res = TRUE` (default) and `save_path` to a string
#' that specifies both the location and the beginning of the file name 
#' (e.g., "~/Documents/run"). The file name will have "_swolca_adapt.RData" or 
#' "_swolca_results.RData" appended to it.
#'
#' If hyperparameters for the adaptive or fixed sampler are left as `NULL` 
#' (default), the following default values are used. Let \eqn{K} refer to 
#' `K_max` for the adaptive sampler and `K_fixed` for the fixed sampler. 
#' For \eqn{\pi}, a Dirichlet prior with hyperparameter \eqn{\alpha = 1/K} for 
#' each component. For \eqn{\theta_{jk\cdot}}, a Dirichlet prior with 
#' hyperparameter  \eqn{\eta_j} equal to `rep(1, R_j)` where `R_j` is the number 
#' of categories for exposure item j. If `R_j < R`, the remaining categories have
#' hyperparameter set to 0.01. This is done independently for each exposure item j
#' and is assumed to be the same across latent classes. For \eqn{\xi_{k\cdot}}, a 
#' Multivariate Normal distribution with mean vector hyperparameter \eqn{\mu_0} 
#' drawn from a Normal(0,1) hyperprior for each component, and variance matrix 
#' hyperparameter \eqn{\Sigma_0} a diagonal matrix with diagonal components drawn 
#' from InvGamma(shape=3.5, scale=6.25) distributions. Note that hyperparameters
#' for the fixed sampler should probably only be specified if running the 
#' fixed sampler directly, bypassing the adaptive sampler. 
#' 
#' @return
#' If the fixed sampler is run, returns an object `res` of class `"swolca"`; a 
#' list containing the following:
#' \describe{
#'   \item{\code{estimates_unadj}}{List of unadjusted posterior model results}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{List of data variables used}
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler}
#'   \item{\code{estimates}}{If `adjust_var = TRUE` (default), list of adjusted 
#'   posterior model results with correct uncertainty estimation}
#' }
#'
#' If `save_res = TRUE` (default), also saves `res` as 
#' `[save_path]_swolca_results.RData`. 
#' 
#' If only the adaptive sampler is run (i.e., `run_sampler` = `"adapt"`), returns
#' list `res` containing:
#' \describe{
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler, 
#' obtained from the adaptive sampler}
#'   \item{\code{K_MCMC}}{Adaptive sampler MCMC output for K}
#' }
#' 
#' If `save_res = TRUE` (default), also saves `res` as 
#' `[save_path]_swolca_adapt.RData`. 
#' 
#' @seealso [solca()] [wolca()]
#' 
#' @importFrom RcppTN rtn
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats rnorm median model.matrix as.formula
#' @export
#'
#' @examples
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
#'        cluster_id = cluster_id, stratum_id = stratum_id, V_data = V_data,
#'        run_sampler = "both", glm_form = glm_form, adapt_seed = 1,
#'        n_runs = 50, burn = 25, thin = 1, save_res = FALSE)
#'    
#' \dontrun{        
#' # Run swolca on NHANES data
#' data("data_nhanes")
#' x_mat <- as.matrix(dplyr::select(data_nhanes, citrus:drinks))
#' y_all <- data_nhanes$BP_flag
#' stratum_id <- data_nhanes$stratum_id
#' cluster_id <- data_nhanes$cluster_id
#' sampling_wt <- data_nhanes$sample_wt
#' V_data <- dplyr::select(data_nhanes, age_cat, racethnic, smoker, physactive)
#' glm_form <- "~ age_cat + racethnic + smoker + physactive"
#' res_nhanes <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
#'                      cluster_id = cluster_id, stratum_id = stratum_id,
#'                      V_data = V_data, run_sampler = "both",
#'                      glm_form = glm_form, adapt_seed = 20230225,
#'                      n_runs = 20000, burn = 19800, thin = 5, save_res = FALSE,
#'                      save_path = "~/Documents/run")
#' }
#'
swolca <- function(x_mat, y_all, sampling_wt, cluster_id = NULL, 
                   stratum_id = NULL, V_data = NULL, run_sampler = "both", 
                   glm_form,
                   K_max = 30, adapt_seed = NULL, class_cutoff = 0.05,
                   alpha_adapt = NULL, eta_adapt = NULL,
                   mu0_adapt = NULL, Sig0_adapt = NULL,
                   fixed_seed = NULL, K_fixed = NULL, 
                   alpha_fixed = NULL, eta_fixed = NULL,
                   mu0_fixed = NULL, Sig0_fixed = NULL,
                   n_runs = 20000, burn = 10000, thin = 5, adjust_var = TRUE,
                   num_reps = 100, save_res = TRUE, save_path = NULL) {
  
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
  
  # Obtain normalized weights
  kappa <- sum(sampling_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
  w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
  
  # If no additional covariates, set V_data to be a column of all ones
  if (is.null(V_data)) {
    V_data <- as.data.frame(matrix(1, nrow = n))
  }
  
  #================= Catch errors ==============================================
  catch_errors(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
               cluster_id = cluster_id, stratum_id = stratum_id, V_data = V_data,
               run_sampler = run_sampler, glm_form = glm_form,
               K_max = K_max, class_cutoff = class_cutoff,
               adapt_seed = adapt_seed, fixed_seed = fixed_seed,
               alpha_adapt = alpha_adapt, eta_adapt = eta_adapt, 
               mu0_adapt = mu0_adapt, Sig0_adapt = Sig0_adapt, 
               K_fixed = K_fixed, alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
               mu0_fixed = mu0_fixed, Sig0_fixed = Sig0_fixed,
               n_runs = n_runs, burn = burn, thin = thin, 
               save_res = save_res, save_path = save_path)

  # Obtain probit regression design matrix without class assignment
  V <- model.matrix(as.formula(glm_form), data = V_data)
  # Number of regression covariates excluding class assignment
  q <- ncol(V)  

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
    # Default hyperparameters for xi
    if (is.null(mu0_adapt)) {
      mu0_adapt <- vector("list", K_max)
      for (k in 1:K_max) {
        # MVN(0,1) hyperprior for prior mean of xi
        mu0_adapt[[k]] <- stats::rnorm(n = q)
      }
    }
    if (is.null(Sig0_adapt)) {
      Sig0_adapt <- vector("list", K_max)
      for (k in 1:K_max) {
        # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
        # components and mean variance 2.5 for a weakly informative prior on xi
        Sig0_adapt[[k]] <- diag(LaplacesDemon::rinvgamma(n = q, shape = 3.5, scale = 6.25),
                                nrow = q, ncol = q)
      }
    }
    
    #================= Initialize OLCA model =====================================
    # Obtain pi, theta, c_all
    OLCA_params <- init_OLCA(alpha = alpha_adapt, eta = eta_adapt, n = n,
                             K = K_max, J = J, R = R)
    
    #================= Initialize probit model ===================================
    # Obtain xi, z_all
    probit_params <- init_probit(mu0 = mu0_adapt, Sig0 = Sig0_adapt, K = K_max,
                                 q = q, n = n, V = V, y_all = y_all,
                                 c_all = OLCA_params$c_all)
    
    #================= Run adaptive sampler to obtain number of classes ========
    # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
    MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params,
                              probit_params = probit_params,
                              n_runs = n_runs, burn = burn, thin = thin,
                              K = K_max, J = J, R = R, n = n, q = q,
                              w_all = w_all, x_mat = x_mat, y_all = y_all, V = V,
                              alpha = alpha_adapt, eta = eta_adapt,
                              Sig0 = Sig0_adapt, mu0 = mu0_adapt)
    
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
    rm(OLCA_params, probit_params, MCMC_out)
  }

  #================= FIXED SAMPLER =============================================
  if (run_sampler %in% c("both", "fixed")) {
    print("Running fixed sampler...")
    
    # Catch errors: check hyperparameter dimensions for fixed sampler
    catch_errors(x_mat = x_mat, K_fixed = K_fixed, 
                 alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
                 mu0_fixed = mu0_fixed, Sig0_fixed = Sig0_fixed,
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
      # Hyperparameter for prior for theta
      # Unviable categories have value 0.01 to prevent rank deficiency issues
      eta_fixed <- matrix(0.01, nrow = J, ncol = R) 
      for (j in 1:J) {
        eta_fixed[j, 1:R_j[j]] <- rep(1, R_j[j]) 
      }
    }
    # Default hyperparameters for xi
    if (is.null(mu0_fixed)) {
      mu0_fixed <- vector("list", K_fixed)
      for (k in 1:K_fixed) {
        # MVN(0,1) hyperprior for prior mean of xi
        mu0_fixed[[k]] <- stats::rnorm(n = q)
      }
    }
    if (is.null(Sig0_fixed)) {
      Sig0_fixed <- vector("list", K_fixed)
      for (k in 1:K_fixed) {
        # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
        # components and mean variance 2.5 for a weakly informative prior on xi
        Sig0_fixed[[k]] <- diag(LaplacesDemon::rinvgamma(n = q, shape = 3.5, scale = 6.25),
                                nrow = q, ncol = q)
      }
    }
    
    #================= Run fixed sampler to obtain posteriors ====================
    # Initialize OLCA model using fixed number of classes. Obtain pi, theta, c_all
    OLCA_params <- init_OLCA(alpha = alpha_fixed, eta = eta_fixed, n = n,
                             K = K_fixed, J = J, R = R)
    
    # Initialize probit model using fixed number of classes. Obtain xi, z_all
    probit_params <- init_probit(mu0 = mu0_fixed, Sig0 = Sig0_fixed, K = K_fixed,
                                 q = q, n = n, V = V, y_all = y_all,
                                 c_all = OLCA_params$c_all)
    
    # Run MCMC algorithm using fixed number of classes
    # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
    MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params,
                              probit_params = probit_params,
                              n_runs = n_runs, burn = burn, thin = thin,
                              K = K_fixed, J = J, R = R, n = n, q = q,
                              w_all = w_all, x_mat = x_mat, y_all = y_all, V = V,
                              alpha = alpha_fixed, eta = eta_fixed,
                              Sig0 = Sig0_fixed, mu0 = mu0_fixed)
    
    # Post-processing to recalibrate labels and remove extraneous empty classes
    # Obtain K_med, pi, theta, xi, loglik_MCMC
    post_MCMC_out <- post_process(MCMC_out = MCMC_out, J = J, R = R, q = q,
                                  class_cutoff = class_cutoff)
    
    # Obtain posterior estimates, reduce number of classes, get estimates
    # Obtain K_red, pi_red, theta_red, xi_red, pi_med, theta_med, xi_med, Phi_med,
    # c_all, pred_class_probs, loglik_med
    estimates <- get_estimates(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out,
                               n = n, J = J, V = V, y_all = y_all, x_mat = x_mat)
    
    # Create output list. Replaces adaptive sampler output list
    res <- list(estimates_unadj = estimates, V = V, MCMC_out = MCMC_out,
                post_MCMC_out = post_MCMC_out, K_fixed = K_fixed)
    
    #================= VARIANCE ADJUSTMENT =====================================
    if (adjust_var) {
      print("Running variance adjustment")
      
      # Stan model
      mod_stan <- stanmodels$SWOLCA_main
      
      # Apply variance adjustment for correct coverage
      # Obtain pi_red, theta_red, xi_red, pi_med, theta_med, xi_med, Phi_med, 
      # c_all, pred_class_probs, log_lik_med
      estimates_adj <- var_adjust(mod_stan = mod_stan, estimates = estimates,
                                  K = estimates$K_red, J = J, R_j = R_j, R = R, 
                                  n = n, q = q, x_mat = x_mat, y_all = y_all, 
                                  V = V, w_all = w_all, alpha = alpha_fixed, 
                                  eta = eta_fixed, mu0 = mu0_fixed, 
                                  Sig0 = Sig0_fixed, stratum_id = stratum_id, 
                                  cluster_id = cluster_id, num_reps = num_reps)
      res$estimates <- estimates_adj 
    }
  }  

  #================= SAVE AND RETURN OUTPUT ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  res$runtime <- runtime

  # Store data variables used
  data_vars <- list(n = n, J = J, R = R, q = q, sample_wt = sampling_wt,
                    X_data = x_mat, Y_data = y_all, V_data = V_data, glm_form = glm_form,
                    true_Si = stratum_id, cluster_id = cluster_id)
  res$data_vars <- data_vars
  
  class(res) <- "swolca"
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_swolca_results.RData"))
  }

  return(res)
}


