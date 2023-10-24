#===================================================
## Supervised Overfitted Latent Class Model
## Programmer: SM Wu   
## Data: Simulations       
#===================================================

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
scen_samp <- args[[1]]  # Simulation scenario
iter_pop <- args[[2]]   # Population iteration
samp_n <- args[[3]]     # Sample number
covs <- args[[4]]       # Probit covariates
if (covs == 1) {
  covs <- NULL
} else if (covs == 2) {
  covs <- "true_Si"
} else if (covs == 3) {
  covs <- "additional"
}

# Load libraries
library(R.matlab)
library(stringr)
library(plyr)
library(dplyr)
library(LaplacesDemon)
library(truncnorm)
library(fastDummies)
library(matrixStats)
library(Matrix)
library(gtools)
library(e1071)
library(rstan)
library(survey)
library(Rcpp)
library(RcppArmadillo)
library(RcppTN)

#========================= MAIN FUNCTION =======================================

# 'SOLCA_main_Rcpp' runs the SOLCA model and saves and returns results
# Inputs:
#   data_path: String path for input dataset
#   adapt_path: String path for adaptive sampler file
#   res_path: String path for output file
#   save_res: Boolean specifying if results should be saved. Default = TRUE
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
#   covs: String vector of covariates to include in probit model. Default = NULL 
# Outputs: Saves and returns list `res` containing:
#   analysis: List of posterior model results
#   runtime: Total runtime for model
#   data_vars: Input dataset
#   MCMC_out: List of full MCMC output
#   post_MCMC_out: List of MCMC output after relabeling
#   K_MCMC: Adaptive sampler MCMC output for K
# Also saves list `adapt_MCMC` containing:
#   MCMC_out: List of full MCMC output
#   K_fixed: Number of classes to use for fixed sampler; output from adaptive sampler
#   K_MCMC: Adaptive sampler MCMC output for K
SOLCA_main_Rcpp <- function(data_path, adapt_path, res_path, save_res = TRUE, 
                            n_runs, burn, thin, covs = NULL) {
  start_time <- Sys.time()
  
  #================= Read in data ==============================================
  print("Read in data")
  load(data_path)
  data_vars <- sim_data
  # data_vars <- readMat(data_path)$sim.data
  # names(data_vars) <- str_replace_all(dimnames(data_vars)[[1]], "[.]", "_")
  
  # Obtain dimensions
  n <- dim(data_vars$X_data)[1]        # Number of individuals
  p <- dim(data_vars$X_data)[2]        # Number of exposure items
  d <- max(apply(data_vars$X_data, 2,  # Number of exposure categories 
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
  # Obtain data
  x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
  y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
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
  
  # Set normalized weights to 1
  w_all <- rep(1, n)
  
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
  # # Post-processing to recalibrate labels and remove extraneous empty classes
  # # Obtain K_med, pi, theta, xi, loglik_MCMC
  # post_MCMC_out <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)
  # # Identify unique classes using modal exposure categories
  # # Posterior median estimate for theta across iterations
  # theta_med_temp <- apply(post_MCMC_out$theta, c(2, 3, 4), median)
  # # Posterior modal exposure categories for each exposure item and reduced class
  # theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
  # # Identify unique classes
  # unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
  # # Get number of unique classes for fixed sampler
  # K_fixed <- length(unique_classes) 
  
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
  
  runtime <- Sys.time() - start_time
  
  #================= Save and return output ====================================
  res <- list(analysis = analysis, runtime = runtime, 
              data_vars = data_vars, V = V, MCMC_out = MCMC_out, 
              post_MCMC_out = post_MCMC_out, K_MCMC = adapt_MCMC$K_MCMC)
  if (save_res) {
    save(res, file = res_path)
  }
  return(list(res = res, adapt_MCMC = adapt_MCMC))
}



#===================== RUN MAIN SOLCA FUNCTION =================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/SWOLCA/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/SWOLCA/"
data_dir <- "Data/July6/"
res_dir <- "Results/July6/"
model_dir <- "Model_Code/"
model <- "sOFMM"

# # Testing code
# scen_samp <- 111111
# iter_pop <- 1
# samp_n <- 1
# n_runs <- 100
# burn <- 50
# thin <- 5
# covs <- "true_Si"
# save_res <- FALSE

# Define paths
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".RData")   # Input dataset
# data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
#                     "_samp", samp_n, ".mat")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_scen", scen_samp, 
                     "_samp", samp_n, ".RData")  # Output file
res_path <- paste0(wd, res_dir, model, "_results_scen", scen_samp, 
                   "_samp", samp_n, ".RData")  # Output file

# Check if results already exist
already_done <- file.exists(res_path)
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
  print(paste0("Running SOLCA_main for scenario ", scen_samp, ' iter ', 
               iter_pop,' samp ', samp_n))
  results <- SOLCA_main_Rcpp(data_path = data_path, adapt_path = adapt_path,
                             res_path = res_path,
                             save_res = save_res, n_runs = n_runs, 
                             burn = burn, thin = thin, covs = covs)
  print(paste0("Runtime: ", results$res$runtime))
}

