#' Run the SWOLCA model
#'
#' @description
#' `swolca` runs a supervised weighted overfitted latent class analysis (SWOLCA)
#' described in Wu et al. (2023) and saves and returns the results. For proper 
#' variance estimation, please run [swolca_var_adjust()] after.
#'
#' @param x_mat Matrix of multivariate categorical exposures. nxJ
#' @param y_all Vector of binary outcomes. nx1
#' @param sampling_wt Vector of survey sampling weights. nx1. Default is `NULL`, 
#' indicating no sampling weights and setting all weights to 1. 
#' @param cluster_id Vector of individual cluster IDs. nx1. Default is `NULL`,
#' indicating each individual is their own cluster.
#' @param stratum_id Vector of individual stratum IDs. nx1. Default is `NULL`,
#' indicating no stratification.
#' @param V_data Dataframe of additional regression covariates. nxQ. Factor 
#' covariates must be converted to factors. If `NULL` (default), no additional 
#' covariates are to be included. All variables in `glm_form` must 
#' be found in `V_data`.
#' @param glm_form String specifying formula for probit regression, excluding 
#' outcome and latent class. For example, `"~ 1"` for the model with only 
#' latent class as covariates. All variables in `glm_form` must be found in `V_data`.
#' Do not specify interaction terms for latent class by additional covariates, 
#' as these terms are already included. 
#' @param run_sampler String specifying which sampler(s) should be run. Must be 
#' one of `"both"` (default), `"fixed"`, or `"adapt"`. See Details.
#' @param K_max Upper limit for number of classes. Default is 30.
#' @param adapt_seed Numeric seed for adaptive sampler. Default is `NULL`.
#' @param K_fixed True number of classes, if known. Default is `NULL`, as it is
#' not necessary for the adaptive sampler. If bypassing the adaptive sampler and
#' running the fixed sampler directly, need to specify a value here. See Details.
#' @param fixed_seed Numeric seed for fixed sampler. Default is `NULL`.
#' @param class_cutoff Minimum class size proportion when determining number of
#' classes. Default is 0.05.
#' @param n_runs Number of MCMC iterations. Default is 20000.
#' @param burn Number of MCMC iterations to drop as a burn-in period. Default is 10000.
#' @param thin Thinning factor for MCMC iterations. Default is 5.
#' @param update Number specifying that MCMC progress updates should be printed 
#' every `update` iterations. Default is 10.
#' @param save_res Boolean specifying if results should be saved. Default = `TRUE`.
#' @param save_path String specifying directory and file name to save results, 
#' e.g., "~/Documents/run". Default is `NULL`.
#' @param alpha_adapt Adaptive sampler hyperparameter for prior for class 
#' membership probabilities \eqn{\pi}. Default is `NULL` and default values are 
#' used (see Details). If specified, must be (`K_max`)x1. 
#' @param eta_adapt Adaptive sampler hyperparameter for prior for item 
#' level probabilities \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class 
#' \eqn{k}, assumed to be the same across classes. Default is `NULL` and default 
#' values are used (see Details). If specified, must be JxR, where J is the 
#' number of exposure items and R is the maximum number of categories for the exposure.
#' @param mu0_adapt Adaptive sampler mean hyperparameters for regression 
#' coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}. Default is `NULL` and 
#' default values are used (see Details). If specified, must be a list of `K_max` 
#' vectors of dimension Qx1, where Q is the number of regression covariates 
#' excluding latent class assignment.  
#' @param Sig0_adapt Adaptive sampler variance hyperparameters for regression 
#' coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}. Default is `NULL` and 
#' default values are used (see Details). If specified, must be a list of `K_max` 
#' QxQ matrices, where Q is the number of regression covariates excluding latent 
#' class assignment. 
#' @param alpha_fixed Fixed sampler hyperparameter for prior for \eqn{\pi}. Default is 
#' `NULL` and default values are used (see Details). If specified, must be (`K_fixed`)x1. 
#' @param eta_fixed Fixed sampler hyperparameter for prior for item level 
#' probabilities \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class \eqn{k}, 
#' assumed to be the same across classes. Default is `NULL` and default values 
#' are used (see Details). If specified, must be JxR, where J is the number of 
#' exposure items and R is the maximum number of categories for the exposure.
#' @param mu0_fixed Fixed sampler mean hyperparameters for regression coefficients
#' \eqn{\xi_{k\cdot}} for each class \eqn{k}. Default is `NULL` and default 
#' values are used (see Details). If specified, must be a list of `K_fixed` 
#' vectors of dimension Qx1, where Q is the number of regression covariates 
#' excluding latent class assignment. 
#' @param Sig0_fixed Fixed sampler variance hyperparameters for regression 
#' coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}. Default is `NULL` and 
#' default values are used (see Details). If specified, must be a list of 
#' `K_fixed` QxQ matrices, where Q is the number of regression covariates 
#' excluding latent class assignment. 
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
#' To run the post-processing variance adjustment to prevent underestimation of 
#' posterior intervals, run the [swolca_var_adjust()] function after running 
#' [swolca()] (see Examples). 
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
#' of levels for exposure item j. If `R_j < R`, the remaining levels have
#' hyperparameter set to 0.01. This is done independently for each exposure item j
#' and is assumed to be the same across latent classes. For \eqn{\xi_{k\cdot}}, a 
#' Multivariate Normal distribution with mean vector hyperparameter \eqn{\mu_0} 
#' drawn from a Normal(0,1) hyperprior for each component, and variance matrix 
#' hyperparameter \eqn{\Sigma_0} a diagonal matrix with diagonal components drawn 
#' from InvGamma(shape=3.5, scale=6.25) distributions. Note that hyperparameters
#' for the fixed sampler should probably only be specified if running the 
#' fixed sampler directly, bypassing the adaptive sampler.
#'
#' To prevent underflow issues, all \eqn{\theta} and \eqn{\xi} parameters are 
#' restricted to have a minimum value of \eqn{1e-8}. If continuous variables are 
#' incorporated into the binary outcome model, any such variables with standard 
#' deviation greater than 5 may result in errors in the variance adjustment. 
#' To avoid these errors, consider standardizing the variable to have mean 0 and 
#' standard deviation 1 or converting the variable into a categorical form. 
#' 
#' @return
#' If the fixed sampler is run, returns an object `res` of class `"swolca"`; a 
#' list containing the following:
#' \describe{
#'   \item{\code{estimates}}{List of posterior model results, resulting from a 
#'   call to [get_estimates()]}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{List of data variables used, including:
#'   `n`: Sample size.
#'   `J`: Number of exposure items.
#'   `R_j`: Vector of number of exposure levels for each item; Jx1.
#'   `R`: Maximum number of exposure categories across items.
#'   `Q`: Number of regression covariates excluding class assignment.
#'   `w_all`: Vector of sampling weights normalized to sum to n; nx1.
#'   `sampling_wt`: Vector of survey sampling weights; nx1.
#'   `x_mat`: Matrix of multivariate categorical exposures; nxJ.
#'   `y_all`: Vector of binary outcomes; nx1.
#'   `V_data`: Dataframe of additional regression covariates; nxQ or NULL. 
#'   `V`: Regression design matrix without class assignment; nxQ.
#'   `glm_form`: String specifying formula for probit regression, excluding 
#' outcome and latent class.
#'   `stratum_id`: Vector of individual stratum IDs; nx1 or NULL.
#'   `cluster_id`: Vector of individual cluster IDs; nx1 or NULL. 
#'   }
#'   \item{\code{MCMC_out}}{List of full MCMC output, resulting from a call to 
#'   [run_MCMC_Rcpp()]}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling, resulting 
#'   from a call to [post_process()]}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler}
#' }
#'
#' If `save_res = TRUE` (default), also saves `res` as 
#' `[save_path]_swolca_results.RData`. 
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
#' 
#' If `save_res = TRUE` (default), also saves `res` as 
#' `[save_path]_swolca_adapt.RData`. 
#' 
#' @seealso [swolca_var_adjust()] [wolca()]
#' 
#' @importFrom RcppTN rtn
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats rnorm median model.matrix as.formula sd
#' @export
#' 
#' @references Wu, S. M., Williams, M. R., Savitsky, T. D., & Stephenson, B. J. 
#' (2023). Derivation of outcome-dependent dietary patterns for low-income women 
#' obtained from survey data using a Supervised Weighted Overfitted Latent Class 
#' Analysis. arXiv preprint arXiv:2310.01575.
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
#'                     cluster_id = cluster_id, stratum_id = stratum_id, V_data = V_data,
#'                     run_sampler = "both", glm_form = glm_form, adapt_seed = 1,
#'                     n_runs = 50, burn = 25, thin = 1, save_res = FALSE)
#'        
#' # Run variance adjustment to prevent underestimation of posterior intervals
#' res_adjust <- swolca_var_adjust(res = res, num_reps = 100, save_res = FALSE, 
#'                                 adjust_seed = 1)    
#'
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
#' res_nhanes_adjust <- swolca_var_adjust(res = res_nhanes, num_reps = 100, 
#'                                        save_res = FALSE, adjust_seed = 1) 
#' }
#'
swolca <- function(x_mat, y_all, sampling_wt = NULL, cluster_id = NULL, 
                   stratum_id = NULL, V_data = NULL, glm_form, 
                   run_sampler = "both", K_max = 30, adapt_seed = NULL, 
                   K_fixed = NULL, fixed_seed = NULL, class_cutoff = 0.05,
                   n_runs = 20000, burn = 10000, thin = 5, update = 10,
                   save_res = TRUE, save_path = NULL,
                   alpha_adapt = NULL, eta_adapt = NULL,
                   mu0_adapt = NULL, Sig0_adapt = NULL,
                   alpha_fixed = NULL, eta_fixed = NULL,
                   mu0_fixed = NULL, Sig0_fixed = NULL) {
  
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
    # Weights norm. constant. If sum(weights)=N, this is 1/(sampling fraction)
    kappa <- sum(sampling_wt) / n  
    # Weights normalized to sum to n, nx1
    w_all <- c(sampling_wt / kappa) 
  }
  
  # If no additional covariates, set V_data to be a column of all ones
  if (is.null(V_data)) {
    V_data <- as.data.frame(matrix(1, nrow = n))
  }
  # Obtain probit regression design matrix without class assignment
  V <- model.matrix(as.formula(glm_form), data = V_data)
  # Number of regression covariates excluding class assignment
  Q <- ncol(V)  
  # Check that continuous covariates do not have variance that is too large
  for (i in 1:ncol(V)) {
    col <- V[, i]
    if(!is.factor(col) & stats::sd(col) > 5) {
      warning(paste0("Standard deviation for continuous covariate ", colnames(V)[i], 
                     " is greater than 5 and may result in estimation errors. ",
                     "To avoid these errors, consider standardizing the variable to have mean 0 ",
                     "and standard deviaiton 1 or converting the variable into a categorical form."))
    }
  }
  
  #================= Catch errors ==============================================
  catch_errors(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
               cluster_id = cluster_id, stratum_id = stratum_id, V_data = V_data,
               run_sampler = run_sampler, glm_form = glm_form, Q = Q,
               K_max = K_max, class_cutoff = class_cutoff,
               adapt_seed = adapt_seed, fixed_seed = fixed_seed,
               alpha_adapt = alpha_adapt, eta_adapt = eta_adapt, 
               mu0_adapt = mu0_adapt, Sig0_adapt = Sig0_adapt, 
               K_fixed = K_fixed, alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
               mu0_fixed = mu0_fixed, Sig0_fixed = Sig0_fixed,
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
    # Default hyperparameters for xi
    if (is.null(mu0_adapt)) {
      mu0_adapt <- vector("list", K_max)
      for (k in 1:K_max) {
        # MVN(0,1) hyperprior for prior mean of xi
        mu0_adapt[[k]] <- stats::rnorm(n = Q)
      }
    }
    if (is.null(Sig0_adapt)) {
      Sig0_adapt <- vector("list", K_max)
      for (k in 1:K_max) {
        # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
        # components and mean variance 2.5 for a weakly informative prior on xi
        Sig0_adapt[[k]] <- diag(LaplacesDemon::rinvgamma(n = Q, shape = 3.5, scale = 6.25),
                                nrow = Q, ncol = Q)
      }
    }
    
    #================= Initialize OLCA model =====================================
    # Obtain pi, theta, c_all
    OLCA_params <- init_OLCA(alpha = alpha_adapt, eta = eta_adapt, n = n,
                             K = K_max, J = J, R = R)
    
    #================= Initialize probit model ===================================
    # Obtain xi, z_all
    probit_params <- init_probit(mu0 = mu0_adapt, Sig0 = Sig0_adapt, K = K_max,
                                 Q = Q, n = n, V = V, y_all = y_all,
                                 c_all = OLCA_params$c_all)
    
    #================= Run adaptive sampler to obtain number of classes ========
    # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
    MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params,
                              probit_params = probit_params,
                              n_runs = n_runs, burn = burn, thin = thin,
                              K = K_max, J = J, R = R, n = n, Q = Q,
                              w_all = w_all, x_mat = x_mat, y_all = y_all, V = V,
                              alpha = alpha_adapt, eta = eta_adapt,
                              Sig0 = Sig0_adapt, mu0 = mu0_adapt, update = update)
    
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
    catch_errors(x_mat = x_mat, K_fixed = K_fixed, Q = Q,
                 alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
                 mu0_fixed = mu0_fixed, Sig0_fixed = Sig0_fixed,
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
    # Default hyperparameters for xi
    if (is.null(mu0_fixed)) {
      mu0_fixed <- vector("list", K_fixed)
      for (k in 1:K_fixed) {
        # MVN(0,1) hyperprior for prior mean of xi
        mu0_fixed[[k]] <- stats::rnorm(n = Q)
      }
    }
    if (is.null(Sig0_fixed)) {
      Sig0_fixed <- vector("list", K_fixed)
      for (k in 1:K_fixed) {
        # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
        # components and mean variance 2.5 for a weakly informative prior on xi
        Sig0_fixed[[k]] <- diag(LaplacesDemon::rinvgamma(n = Q, shape = 3.5, scale = 6.25),
                                nrow = Q, ncol = Q)
      }
    }
    
    #================= Run fixed sampler to obtain posteriors ====================
    # Initialize OLCA model using fixed number of classes. Obtain pi, theta, c_all
    OLCA_params <- init_OLCA(alpha = alpha_fixed, eta = eta_fixed, n = n,
                             K = K_fixed, J = J, R = R)
    
    # Initialize probit model using fixed number of classes. Obtain xi, z_all
    probit_params <- init_probit(mu0 = mu0_fixed, Sig0 = Sig0_fixed, K = K_fixed,
                                 Q = Q, n = n, V = V, y_all = y_all,
                                 c_all = OLCA_params$c_all)
    
    # Run MCMC algorithm using fixed number of classes
    # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
    MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params,
                              probit_params = probit_params,
                              n_runs = n_runs, burn = burn, thin = thin,
                              K = K_fixed, J = J, R = R, n = n, Q = Q,
                              w_all = w_all, x_mat = x_mat, y_all = y_all, V = V,
                              alpha = alpha_fixed, eta = eta_fixed,
                              Sig0 = Sig0_fixed, mu0 = mu0_fixed, update = update)
    
    # Post-processing to recalibrate labels and remove extraneous empty classes
    # Obtain K_med, pi, theta, xi, loglik_MCMC
    post_MCMC_out <- post_process(MCMC_out = MCMC_out, J = J, R = R, Q = Q,
                                  class_cutoff = class_cutoff)
    
    # Obtain posterior estimates, reduce number of classes, get estimates
    # Obtain K_red, pi_red, theta_red, xi_red, pi_med, theta_med, xi_med, Phi_med,
    # c_all, pred_class_probs, loglik_med
    estimates <- get_estimates(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out,
                               n = n, J = J, V = V, y_all = y_all, x_mat = x_mat)
    
    # Create output list. Replaces adaptive sampler output list
    res <- list(estimates = estimates, MCMC_out = MCMC_out,
                post_MCMC_out = post_MCMC_out, K_fixed = K_fixed)
  }  

  #================= SAVE AND RETURN OUTPUT ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  res$runtime <- runtime

  # Store data variables used
  data_vars <- list(n = n, J = J, R_j = R_j, R = R, Q = Q, w_all = w_all, 
                    sampling_wt = sampling_wt, x_mat = x_mat, y_all = y_all, 
                    V_data = V_data, V = V, glm_form = glm_form, 
                    stratum_id = stratum_id, cluster_id = cluster_id)
  res$data_vars <- data_vars
  
  class(res) <- "swolca"
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_swolca_results.RData"))
  }

  return(res)
}


