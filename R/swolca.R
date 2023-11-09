#' Run the SWOLCA model
#'
#' @description
#' `swolca` runs a supervised weighted overfitted latent class analysis (SWOLCA)
#' and saves and returns the results.
#'
#' @param x_mat Categorical exposure matrix. nxp
#' @param y_all Vector of outcomes. nx1
#' @param sampling_wt Vector of survey sampling weights. nx1. If no sampling 
#' weights are available, set this to a vector of ones. 
#' @param cluster_id Vector of individual cluster IDs. nx1. Default is `NULL`,
#' indicating each individual is their own cluster.
#' @param stratum_id Vector of individual stratum IDs. nx1. Default is `NULL`,
#' indicating no stratification.
#' @param V Regression design matrix without class assignment. nxq
#' @param run_sampler String specifying which sampler(s) should be run. Must be 
#' one of `"both"` (default), `"fixed"`, or `"adapt"`.
#' @param K_max Upper limit for number of classes. Default is 30.
#' @param adapt_seed Numeric seed for adaptive sampler. Default is `NULL`.
#' @param class_cutoff Minimum class size proportion when determining number of
#' classes in adaptive sampler. Default is 0.05.
#' @param alpha_adapt Adaptive sampler hyperparameter for prior for class 
#' membership probabilities \eqn{\pi}. Default is `NULL` and default values are 
#' used. If specified, must be (`K_max`)x1. 
#' @param eta_adapt Adaptive sampler hyperparameter for prior for item consumption
#' level probabilities \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class 
#' \eqn{k}. Default is `NULL` and default values are used. If specified, must be 
#' Rx1, where R is the maximum number of categories for the exposure.
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
#' \eqn{k}. Default is `NULL` and default values are used. If specified, must be 
#' Rx1, where R is the maximum number of categories for the exposure.
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
#' @param save_res Boolean specifying if results should be saved. Default = `TRUE`.
#' @param save_path String specifying directory and file name to save results. Default is `NULL`.
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
#' (e.g., "~/Documents/run"). The file name will have "_swolca_adapt.RData" or 
#' "_swolca_results.RData" appended to it.
#'
#' If hyperparameters for the adaptive or fixed sampler are left as `NULL` 
#' (default), the following default values are used. Let \eqn{K} refer to 
#' `K_max` for the adaptive sampler and `K_fixed` for the fixed sampler. 
#' For \eqn{\pi}, a Dirichlet prior with hyperparameter \eqn{\alpha = 1/K} for 
#' each component. For \eqn{\theta_{jk\cdot}}, a Dirichlet prior with 
#' hyperparameter \eqn{\eta = 1} for each component. For \eqn{\xi_{k\cdot}}, a 
#' Multivariate Normal distribution with mean vector hyperparameter \eqn{\mu_0} 
#' drawn from a Normal(0,1) hyperprior for each component, and variance matrix 
#' hyperparameter \eqn{\Sigma_0} a diagonal matrix with diagonal components drawn 
#' from InvGamma(shape=5/2, scale=2/5) distributions. Note that hyperparameters
#' for the fixed sampler should probably only be specified if running the 
#' fixed sampler directly, bypassing the adaptive sampler. 
#' 
#' @return
#' If the fixed sampler is run, returns list `res` containing:
#' \describe{
#'   \item{\code{estimates}}{List of adjusted posterior model results}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{List of data variables used}
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler}
#'   \item{\code{estimates_unadj}}{List of unadjusted posterior model results}
#'   \item{\code{K_MCMC}}{If `K_fixed = NULL` and the adaptive sampler is run,
#'   output list also contains MCMC output for the number of classes with size
#'   greater than `class_cutoff` for each iteration}
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
#' @importFrom stats rnorm median
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
#' V <- matrix(1, nrow = n) # Regression design matrix without class assignment
#' 
#' # Run swolca
#' res <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
#'        cluster_id = cluster_id, stratum_id = stratum_id, V = V, 
#'        run_sampler = "both", adapt_seed = 1,
#'        n_runs = 50, burn = 25, thin = 1, save_res = FALSE)
#'
swolca <- function(x_mat, y_all, sampling_wt, cluster_id = NULL, 
                   stratum_id = NULL, V, run_sampler = "both",
                   K_max = 30, adapt_seed = NULL, class_cutoff = 0.05,
                   alpha_adapt = NULL, eta_adapt = NULL,
                   mu0_adapt = NULL, Sig0_adapt = NULL,
                   fixed_seed = NULL, K_fixed = NULL, 
                   alpha_fixed = NULL, eta_fixed = NULL,
                   mu0_fixed = NULL, Sig0_fixed = NULL,
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
  q <- ncol(V)              # Number of regression covariates excluding class assignment
  
  # Obtain normalized weights
  kappa <- sum(sampling_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
  w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
  
  #================= Catch errors ==============================================
  catch_errors(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
               cluster_id = cluster_id, stratum_id = stratum_id, V = V,
               run_sampler = run_sampler, 
               K_max = K_max, class_cutoff = class_cutoff,
               alpha_adapt = alpha_adapt, eta_adapt = eta_adapt, 
               mu0_adapt = mu0_adapt, Sig0_adapt = Sig0_adapt, 
               K_fixed = K_fixed, alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
               mu0_fixed = mu0_fixed, Sig0_fixed = Sig0_fixed,
               n_runs = n_runs, burn = burn, thin = thin, 
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
      eta_adapt <- rep(1, R)                 # Hyperparameter for prior for theta
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
      eta_fixed <- rep(1, R)                 # Hyperparameter for prior for theta
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
    
    
    #================= VARIANCE ADJUSTMENT =======================================
    print("Running variance adjustment")
    
    # Stan model
    mod_stan <- stanmodels$SWOLCA_main
    
    # Apply variance adjustment for correct coverage
    # Obtain pi_red_adj, theta_red_adj, xi_red_adj, pi_med_adj, theta_med_adj,
    # xi_med_adj, Phi_med_adj, c_all, pred_class_probs, log_lik_med
    estimates_adj <- var_adjust(mod_stan = mod_stan, estimates = estimates,
                                K = estimates$K_red, J = J, R = R, n = n, q = q,
                                x_mat = x_mat, y_all = y_all, V = V, w_all = w_all,
                                stratum_id = stratum_id, cluster_id = cluster_id)
    
    # Create output list. Replaces adaptive sampler output list
    res <- list(estimates = estimates_adj, V = V, MCMC_out = MCMC_out,
                post_MCMC_out = post_MCMC_out, K_fixed = K_fixed,
                estimates_unadj = estimates)
  }  

  #================= SAVE AND RETURN OUTPUT ====================================
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
    save(res, file = paste0(save_path, "_swolca_results.RData"))
  }

  return(res)
}


