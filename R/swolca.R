#' Run the SWOLCA model
#'
#' @description
#' `swolca` runs a supervised weighted overfitted latent class estimates (SWOLCA)
#' and saves and returns the results.
#'
#' @param x_mat Categorical exposure matrix. nxp
#' @param y_all Vector of outcomes. nx1
#' @param sampling_wt Vector of survey sampling weights. nx1
#' @param cluster_id Vector of individual cluster IDs. nx1
#' @param stratum_id Vector of individual stratum IDs. nx1
#' @param V Regression design matrix without class assignment. nxq
#' @param K_max Upper limit for number of classes. Default is 30.
#' @param adapt_seed Numeric seed for adaptive sampler. Default is `NULL`.
#' @param class_cutoff Minimum class size proportion when determining number of
#' classes in adaptive sampler. Default is 0.05.
#' @param alpha_adapt Adaptive sampler hyperparameter for prior for pi. Default is `NULL`.
#' @param eta_adapt Adaptive sampler hyperparameter for prior for theta. Default is `NULL`.
#' @param mu0_adapt List of K qx1 vectors of adaptive sampler mean hyperparameters for xi. Default is `NULL`.
#' @param Sig0_adapt List of K qxq matrices of adaptive sampler variance hyperparameters for xi. Default is `NULL`.
#' @param alpha_fixed Fixed sampler hyperparameter for prior for pi. Default is `NULL`.
#' @param eta_fixed Fixed sampler hyperparameter for prior for theta. Default is `NULL`.
#' @param mu0_fixed List of K qx1 vectors of fixed sampler mean hyperparameters for xi. Default is `NULL`.
#' @param Sig0_fixed List of K qxq matrices of fixed sampler variance hyperparameters for xi. Default is `NULL`.
#' @param K_true True number of classes, if known. If `NULL` (default),
#' adaptive sampler is run. Otherwise, model uses fixed sampler directly.
#' @param fixed_seed Numeric seed for fixed sampler. Default is `NULL`.
#' @param n_runs Number of MCMC iterations. Default is 20000.
#' @param burn Number of MCMC iterations to drop as a burn-in period. Default is 10000.
#' @param thin Thinning factor for MCMC iterations. Default is 5.
#' @param mod_stan Stan model to be used for variance adjustment
#' @param save_res Boolean specifying if results should be saved. Default = `TRUE`.
#' @param save_path String specifying directory and file name to save results. Default is `NULL`.
#' 
#' @details 
#' If hyperparameters are left as `NULL` (default), the following values are 
#' used. 
#'
#' @return
#' Saves and returns list `res` containing:
#' \describe{
#'   \item{\code{estimates_adj}}{List of adjusted posterior model results}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{List of data variables used}
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling}
#'   \item{\code{K_MCMC}}{Adaptive sampler MCMC output for K}
#' }
#'
#' Also saves list `adapt_MCMC` containing:
#' \describe{
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{K_fixed}}{Number of classes to use for fixed sampler; output
#'   from adaptive sampler}
#'   \item{\code{K_MCMC}}{Adaptive sampler MCMC output for K}
#' }
#'
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' cluster_id <- data_vars$cluster_id  # Cluster indicators, nx1
#' sampling_wt <- data_vars$sample_wt
#' # Probit model includes C and S: C + S + C:S
#' stratum_id <- data_vars[[covs]]  # Stratifying variable, nx1
#' # Regression design matrix without class assignment, nxq
#' V_data <- data.frame(s = as.factor(stratum_id))
#' V <- model.matrix(~ s, V_data)
#' q <- ncol(V)  # Number of regression covariates excluding class assignment
#'
#' # Stan model
#' mod_stan <- stanmodels$WSOLCA_main
#' # Run swolca
#' swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
#'        cluster_id = cluster_id, stratum_id = stratum_id, V = V, adapt_seed = 1,
#'        n_runs = 200, burn = 100, thin = 5, mod_stan = mod_stan, save_res = FALSE)
#'
swolca <- function(x_mat, y_all, sampling_wt, cluster_id, stratum_id, V,
                   K_max = 30, adapt_seed = NULL, class_cutoff = 0.05,
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
  n <- dim(x_mat)[1]        # Number of individuals
  p <- dim(x_mat)[2]        # Number of exposure items
  d <- max(apply(x_mat, 2,  # Number of exposure categories
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
  q <- ncol(V)

  # Obtain normalized weights
  kappa <- sum(sampling_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
  w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1

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
        mu0_adapt[[k]] <- rnorm(n = q)
      }
    }
    if (is.null(Sig0_adapt)) {
      Sig0_adapt <- vector("list", K_max)
      for (k in 1:K_max) {
        # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
        # components and mean variance 2.5 for a weakly informative prior on xi
        Sig0_adapt[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25),
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
    K_med <- round(median(K_MCMC))
    # Get number of unique classes for fixed sampler
    K_fixed <- K_med
    print(paste0("K_fixed: ", K_fixed))
    # Save adaptive output
    adapt_MCMC <- list(MCMC_out = MCMC_out, K_fixed = K_fixed, K_MCMC = K_MCMC)
    if (save_res) {
      save(adapt_MCMC, file = paste0(save_path, "adapt.RData"))
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
      mu0_fixed[[k]] <- rnorm(n = q)
    }
  }
  if (is.null(Sig0_fixed)) {
    Sig0_fixed <- vector("list", K_fixed)
    for (k in 1:K_fixed) {
      # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
      # components and mean variance 2.5 for a weakly informative prior on xi
      Sig0_fixed[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25),
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


  #================= VARIANCE ADJUSTMENT =======================================
  print("Variance adjustment")
  # Create Stan model

  # Apply variance adjustment for correct coverage
  # Obtain pi_red_adj, theta_red_adj, xi_red_adj, pi_med_adj, theta_med_adj,
  # xi_med_adj, Phi_med_adj, c_all, pred_class_probs, log_lik_med
  estimates_adj <- var_adjust(mod_stan = mod_stan, estimates = estimates,
                             K = estimates$K_red, p = p, d = d, n = n, q = q,
                             x_mat = x_mat, y_all = y_all, V = V, w_all = w_all,
                             stratum_id = stratum_id, cluster_id = cluster_id)


  #================= SAVE AND RETURN OUTPUT ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time

  # Store data variables used
  data_vars <- list(n = n, p = p, d = d, q = q, sample_wt = sampling_wt,
                    X_data = x_mat, Y_data = y_all, V = V,
                    true_Si = stratum_id, cluster_id = cluster_id)
  # Create output list
  res <- list(estimates_adj = estimates_adj, runtime = runtime,
              data_vars = data_vars, V = V, MCMC_out = MCMC_out,
              post_MCMC_out = post_MCMC_out, K_MCMC = adapt_MCMC$K_MCMC)
  if (save_res) {
    save(res, file = paste0(save_path, "results.RData"))
  }

  return(list(res = res, adapt_MCMC = adapt_MCMC))
}


