#' Run the SOLCA model
#'
#' @description
#' `solca` runs an unweighted supervised overfitted latent class analysis (SOLCA)
#' and saves and returns the results.
#' 
#' @inheritParams swolca
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
#' If `save_res = TRUE` (default), also saves `res` as `[save_path]_solca_results.RData`
#' and, if `K_true = NULL` so that the adaptive sampler is run, list `adapt_MCMC`
#' is saved as  `[save_path]_solca_adapt.RData`. List `adapt_MCMC` contains:
#' \describe{
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler, 
#' obtained from the adaptive sampler}
#'   \item{\code{K_MCMC}}{Adaptive sampler MCMC output for K}
#' }
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
#' n <- dim(x_mat)[1]        # Number of individuals
#' 
#' # Probit model only includes latent class
#' V <- matrix(1, nrow = n) # Regression design matrix without class assignment
#' 
#' # Run solca
#' res <- solca(x_mat = x_mat, y_all = y_all, V = V, adapt_seed = 1, n_runs = 50, 
#'              burn = 25, thin = 1, save_res = FALSE)
#'
solca <- function(x_mat, y_all, V, K_max = 30, adapt_seed = NULL, 
                  class_cutoff = 0.05, alpha_adapt = NULL, eta_adapt = NULL,
                  mu0_adapt = NULL, Sig0_adapt = NULL,
                  alpha_fixed = NULL, eta_fixed = NULL,
                  mu0_fixed = NULL, Sig0_fixed = NULL,
                  K_true = NULL, fixed_seed = NULL,
                  n_runs = 20000, burn = 10000, thin = 5, 
                  save_res = TRUE, save_path = NULL) {
  
  # Begin runtime tracker
  start_time <- Sys.time()

  #================= Read in data ==============================================
  print("Read in data")

  # Obtain dimensions
  n <- dim(x_mat)[1]        # Number of individuals
  p <- dim(x_mat)[2]        # Number of exposure items
  d <- max(apply(x_mat, 2,  # Number of exposure categories
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
  q <- ncol(V)              # Number of regression covariates excluding class assignment

  # Set normalized weights to 1
  w_all <- rep(1, n)

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
                             K = K_max, p = p, d = d)
    
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
                              K = K_max, p = p, d = d, n = n, q = q,
                              w_all = w_all, x_mat = x_mat, y_all = y_all, V = V,
                              alpha = alpha_adapt, eta = eta_adapt,
                              Sig0 = Sig0_adapt, mu0 = mu0_adapt)
    
    #================= Post-processing for adaptive sampler ====================
    # Get median number of classes with >= 5% of individuals, over all iterations
    M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
    K_MCMC <- rowSums(MCMC_out$pi_MCMC >= class_cutoff)
    K_med <- round(stats::median(K_MCMC))
    # Get number of unique classes for fixed sampler
    K_fixed <- K_med
    print(paste0("K_fixed: ", K_fixed))
    # Save adaptive output
    adapt_MCMC <- list(MCMC_out = MCMC_out, K_fixed = K_fixed, K_MCMC = K_MCMC)
    if (save_res) {
      save(adapt_MCMC, file = paste0(save_path, "_solca_adapt.RData"))
    }
    # Reduce memory burden
    rm(OLCA_params, probit_params, MCMC_out)
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
                           K = K_fixed, p = p, d = d)
  
  # Initialize probit model using fixed number of classes. Obtain xi, z_all
  probit_params <- init_probit(mu0 = mu0_fixed, Sig0 = Sig0_fixed, K = K_fixed,
                               q = q, n = n, V = V, y_all = y_all,
                               c_all = OLCA_params$c_all)
  
  # Run MCMC algorithm using fixed number of classes
  # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
  MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params,
                            probit_params = probit_params,
                            n_runs = n_runs, burn = burn, thin = thin,
                            K = K_fixed, p = p, d = d, n = n, q = q,
                            w_all = w_all, x_mat = x_mat, y_all = y_all, V = V,
                            alpha = alpha_fixed, eta = eta_fixed,
                            Sig0 = Sig0_fixed, mu0 = mu0_fixed)

  # Post-processing to recalibrate labels and remove extraneous empty classes
  # Obtain K_med, pi, theta, xi, loglik_MCMC
  post_MCMC_out <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)
  
  # Obtain posterior estimates, reduce number of classes, get estimates
  # Obtain K_red, pi_red, theta_red, xi_red, pi_med, theta_med, xi_med, Phi_med,
  # c_all, pred_class_probs, loglik_med
  estimates <- get_estimates(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out,
                             n = n, p = p, V = V, y_all = y_all, x_mat = x_mat)

  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  
  # Store data variables used
  data_vars <- list(n = n, p = p, d = d, q = q, X_data = x_mat, Y_data = y_all, 
                    V = V)
  # Create output list
  res <- list(estimates = estimates, runtime = runtime,
              data_vars = data_vars, V = V, MCMC_out = MCMC_out,
              post_MCMC_out = post_MCMC_out, K_fixed = K_fixed)
  if (is.null(K_true)) {
    res$K_MCMC <- K_MCMC
  }
  if (save_res) {
    save(res, file = paste0(save_path, "_solca_results.RData"))
  }
  
  return(res)
}


