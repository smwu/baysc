#' Run the SWOLCA model
#'
#' @description
#' `swolca` runs a supervised weighted overfitted latent class analysis (SWOLCA)
#' and saves and returns the results.
#'
#' @param data_path String path for input dataset.
#' @param adapt_path String path for adaptive sampler file.
#' @param adj_path String path for adjusted output file.
#' @param stan_path String path for Stan file.
#' @param save_res Boolean specifying if results should be saved. Default = `TRUE`.
#' @param n_runs Number of MCMC iterations.
#' @param burn Whole number specifying burn-in period.
#' @param thin Whole number specifying thinning factor.
#' @param covs String vector of covariates to include in probit model. Default = `NULL`.
#'
#' @return
#' Saves and returns list `res` containing:
#' \describe{
#'   \item{\code{analysis_adj}}{List of adjusted posterior model results}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{Input dataset}
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling}
#'   \item{\code{K_MCMC}}{Adaptive sampler MCMC output for K}
#' }
#'
#' Also saves list `adapt_MCMC` containing:
#' \describe{
#'   \item{\code{MCMC_out}{List of full MCMC output}
#'   \item{\code{K_fixed}{Number of classes to use for fixed sampler; output
#'   from adaptive sampler}
#'   \item{\code{K_MCMC}{Adaptive sampler MCMC output for K}
#'
#' @export
#'
#' @examples
#'

# # Testing code
# scen_samp <- 111211
# iter_pop <- 1
# samp_n <- 1
#
# n_runs <- 500
# burn <- 250
# thin <- 5
# covs <- "true_Si"
# save_res <- FALSE

# Define paths
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".RData")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_comb_scen", scen_samp,
                     "_samp", samp_n, ".RData")  # Output file
adj_path <- paste0(wd, res_dir, model, "_results_comb_adjRcpp_scen", scen_samp,
                   "_samp", samp_n, ".RData")      # Adjusted output file
stan_path <- paste0(wd, model_dir, "WSOLCA_main.stan")  # Stan file

# Check if results already exist
already_done <- file.exists(adj_path)
if (already_done) {
  print(paste0('Scenario ', scen_samp, ' iter ', iter_pop, ' samp ', samp_n,
               ' already exists.'))
} else {
  n_runs <- 20000
  burn <- 10000
  thin <- 5
  save_res <- TRUE
  # covs <- "true_Si"

  # Source R helper functions
  source(paste0(wd, model_dir, "helper_functions.R"))
  # Source Rcpp functions
  Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))
  # Set seed
  set.seed(samp_n)
  # Run model
  print(paste0("Running WSOLCA_main for scenario ", scen_samp, ' iter ',
               iter_pop,' samp ', samp_n))
  results_adj <- WSOLCA_main_Rcpp(data_path = data_path, adapt_path = adapt_path,
                                  adj_path = adj_path, stan_path = stan_path,
                                  save_res = save_res, n_runs = n_runs,
                                  burn = burn, thin = thin, covs = covs)
  print(paste0("Runtime: ", results_adj$res$runtime))
}
WSOLCA_main_Rcpp <- function(data_path, adapt_path, adj_path, stan_path,
                             save_res = TRUE, n_runs, burn, thin,
                             covs = NULL) {
  start_time <- Sys.time()

  #================= Read in data ==============================================
  print("Read in data")
  load(data_path)
  data_vars <- sim_data

  # Obtain dimensions
  n <- dim(data_vars$X_data)[1]        # Number of individuals
  p <- dim(data_vars$X_data)[2]        # Number of exposure items
  d <- max(apply(data_vars$X_data, 2,  # Number of exposure categories
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
  # Obtain data
  x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
  y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
  clus_id_all <- data_vars$cluster_id  # Cluster indicators, nx1
  if (is.null(covs)) {
    # Probit model only includes latent class C
    # No stratifying variable
    s_all <- NULL
    V <- matrix(1, nrow = n)
    q <- 1
  } else if (covs == "true_Si") {
    # Probit model includes C and S: C + S + C:S
    # Stratifying variable, nx1
    s_all <- data_vars[[covs]]
    # Regression design matrix without class assignment, nxq
    V_data <- data.frame(s = as.factor(s_all))
    V <- model.matrix(~ s, V_data)
    # Number of regression covariates excluding class assignment
    q <- ncol(V)
  } else if (covs == "additional") {
    # Probit model includes C, S, A (binary), and B (continuous)
    # C + S + A + B + C:S + C:A + C:B
    # Stratifying variable, nx1
    s_all <- data_vars[["true_Si"]]
    a_all <- data_vars[["true_Ai"]]
    b_all <- data_vars[["true_Bi"]]
    # Regression design matrix without class assignment, nxq
    V_data <- data.frame(s = as.factor(s_all), a = as.factor(a_all), b = b_all)
    V <- model.matrix(~ s + a + b, V_data)
    # Number of regression covariates excluding class assignment
    q <- ncol(V)
  } else {
    stop("Error: covs must be one of 'true_Si', 'additional', or NULL")
  }

  # Obtain normalized weights
  kappa <- sum(data_vars$sample_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
  w_all <- c(data_vars$sample_wt / kappa) # Weights normalized to sum to n, nx1

  #================= ADAPTIVE SAMPLER ==========================================
  print("Adaptive sampler")
  #================= Initialize priors and variables for OLCA model ============
  K_max <- 30                      # Upper limit for number of classes
  alpha <- rep(1, K_max) / K_max   # Hyperparameter for prior for pi
  eta <- rep(1, d)                 # Hyperparameter for prior for theta
  # Obtain pi, theta, c_all
  OLCA_params <- init_OLCA(alpha = alpha, eta = eta, n = n, K = K_max, p = p,
                           d = d)

  #================= Initialize priors and variables for probit model ==========
  # Initialize hyperparameters for xi
  mu0 <- Sig0 <- vector("list", K_max)
  for (k in 1:K_max) {
    # MVN(0,1) hyperprior for prior mean of xi
    mu0[[k]] <- rnorm(n = q)
    # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
    # components and mean variance 2.5 for a weakly informative prior on xi
    Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
  }
  # Obtain xi, z_all
  probit_params <- init_probit(mu0 = mu0, Sig0 = Sig0, K = K_max, q = q, n = n,
                               V = V, y_all = y_all, c_all = OLCA_params$c_all)

  #================= Run adaptive sampler to obtain number of classes ==========
  # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
  MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, probit_params = probit_params,
                            n_runs = n_runs, burn = burn,
                            thin = thin, K = K_max, p = p, d = d, n = n, q = q,
                            w_all = w_all, x_mat = x_mat, y_all = y_all, V = V,
                            alpha = alpha, eta = eta, Sig0 = Sig0, mu0 = mu0)

  #================= Post-processing for adaptive sampler ======================
  # Get median number of classes with >= 5% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_MCMC <- rowSums(MCMC_out$pi_MCMC >= 0.05)
  K_med <- round(median(K_MCMC))
  # Get number of unique classes for fixed sampler
  K_fixed <- K_med
  print(paste0("K_fixed: ", K_fixed))
  # Save adaptive output
  adapt_MCMC <- list(MCMC_out = MCMC_out, K_fixed = K_fixed, K_MCMC = K_MCMC)
  if (save_res) {
    save(adapt_MCMC, file = adapt_path)
  }
  # Reduce memory burden
  rm(OLCA_params, probit_params, MCMC_out)

  # # Obtain K_med, pi, theta, xi
  # post_MCMC_out <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)
  #
  # # Identify unique classes using modal exposure categories
  # # Posterior median estimate for theta across iterations
  # theta_med_temp <- apply(post_MCMC_out$theta, c(2, 3, 4), median)
  # # Posterior modal exposure categories for each exposure item and reduced class
  # theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
  # # Identify unique classes
  # unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
  #
  # # Get number of unique classes for fixed sampler
  # K_fixed <- length(unique_classes)
  #
  # # Reduce memory burden
  # rm(OLCA_params, probit_params, MCMC_out, post_MCMC_out)

  #================= FIXED SAMPLER =============================================
  print("Fixed sampler")
  set.seed(20230629)
  #================= Run fixed sampler to obtain posteriors ====================
  # Initialize OLCA model using fixed number of classes
  # alpha <- rep(2, K_fixed) # Hyperparameter for prior for pi
  alpha <- rep(1, K_fixed) / K_fixed  # Hyperparameter for prior for pi
  # Obtain pi, theta, c_all
  OLCA_params <- init_OLCA(alpha = alpha, eta = eta, n = n, K = K_fixed, p = p,
                           d = d)

  # Initialize probit model using fixed number of classes
  # Initialize hyperparameters for xi
  mu0 <- Sig0 <- vector("list", K_fixed)
  for (k in 1:K_fixed) {
    # MVN(0,1) hyperprior for prior mean of xi
    mu0[[k]] <- rnorm(n = q)
    # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
    # components and mean variance 2.5 for a weakly informative prior on xi
    Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
  }
  # Obtain xi, z_all
  probit_params <- init_probit(mu0 = mu0, Sig0 = Sig0, K = K_fixed, q = q, n = n,
                               V = V, y_all = y_all, c_all = OLCA_params$c_all)

  # Run MCMC algorithm using fixed number of classes
  # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
  MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, probit_params = probit_params,
                            n_runs = n_runs, burn = burn, thin = thin, K = K_fixed,
                            p = p, d = d, n = n, q = q, w_all = w_all, x_mat = x_mat,
                            y_all = y_all, V = V, alpha = alpha, eta = eta,
                            Sig0 = Sig0, mu0 = mu0)

  # Post-processing to recalibrate labels and remove extraneous empty classes
  # Obtain K_med, pi, theta, xi, loglik_MCMC
  post_MCMC_out <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)

  # Obtain posterior estimates, reduce number of classes, analyze results
  # Obtain K_red, pi_red, theta_red, xi_red, pi_med, theta_med, xi_med, Phi_med,
  # c_all, pred_class_probs, loglik_med
  analysis <- analyze_results(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out,
                              n = n, p = p, V = V, y_all = y_all, x_mat = x_mat)

  #================= VARIANCE ADJUSTMENT =======================================
  print("Variance adjustment")
  # Create Stan model
  mod_stan <- stan_model(stan_path)
  # Apply variance adjustment for correct coverage
  # Obtain pi_red_adj, theta_red_adj, xi_red_adj, pi_med_adj, theta_med_adj,
  # xi_med_adj, Phi_med_adj, c_all, pred_class_probs, log_lik_med
  analysis_adj <- var_adjust(mod_stan = mod_stan, analysis = analysis,
                             K = analysis$K_red, p = p, d = d, n = n, q = q,
                             x_mat = x_mat, y_all = y_all, V = V, w_all = w_all,
                             s_all = s_all, clus_id_all = clus_id_all)

  runtime <- Sys.time() - start_time

  #================= Save and return output ====================================
  res <- list(analysis_adj = analysis_adj, runtime = runtime,
              data_vars = data_vars, V = V, MCMC_out = MCMC_out,
              post_MCMC_out = post_MCMC_out, K_MCMC = adapt_MCMC$K_MCMC)
  if (save_res) {
    save(res, file = adj_path)
  }
  return(list(res = res, adapt_MCMC = adapt_MCMC))
}


