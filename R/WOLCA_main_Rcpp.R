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

# 'WOLCA_main_Rcpp' runs the WOLCA model and saves and returns results
# Inputs:
#   data_path: String path for input dataset
#   adapt_path: String path for adaptive sampler file
#   res_path: String path for output file
#   save_res: Boolean specifying if results should be saved. Default = TRUE
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
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
WOLCA_main_Rcpp <- function(data_path, adapt_path, res_path, save_res = TRUE, 
                            n_runs, burn, thin, covs = NULL) {
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
  
  #================= Run adaptive sampler to obtain number of classes ==========
  # Obtain pi_MCMC, theta_MCMC, c_all_MCMC
  MCMC_out <- run_MCMC_Rcpp_WOLCA(OLCA_params = OLCA_params,  
                                  n_runs = round(n_runs/2), burn = round(burn/2), 
                                  thin = thin, K = K_max, p = p, d = d, n = n, 
                                  w_all = w_all, x_mat = x_mat, alpha = alpha, 
                                  eta = eta)
  
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
  rm(OLCA_params, MCMC_out)
  
  
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
  
  # Run MCMC algorithm using fixed number of classes
  # Obtain pi_MCMC, theta_MCMC, c_all_MCMC
  MCMC_out <- run_MCMC_Rcpp_WOLCA(OLCA_params = OLCA_params,  
                                 n_runs = n_runs, burn = burn, thin = thin, K = K_fixed, 
                                 p = p, d = d, n = n, w_all = w_all, x_mat = x_mat, 
                                 alpha = alpha, eta = eta)
  
  # Post-processing to recalibrate labels and remove extraneous empty classes
  # Obtain K_med, pi, theta, xi, loglik_MCMC
  post_MCMC_out <- post_process_WOLCA(MCMC_out = MCMC_out, p = p, d = d)
  
  # Obtain posterior estimates, reduce number of classes, analyze results
  # Obtain K_red, pi_red, theta_red, xi_red, pi_med, theta_med, xi_med, Phi_med, 
  # c_all, pred_class_probs, loglik_med
  analysis <- analyze_results_WOLCA(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out, 
                                  n = n, p = p, x_mat = x_mat)
  
  #================= Fit probit model ==========================================
  if (is.null(covs)) {
    # Probit model only includes latent class C
    svy_data <- data.frame(x = x_mat,
                           y = y_all, 
                           wts = w_all,
                           c = factor(analysis$c_all),
                           clus = clus_id_all)
    svydes <- svydesign(ids = ~clus, weights = ~wts, data = svy_data)
    fit <- svyglm(y ~ c, design = svydes, family = quasibinomial(link = "probit"))
  } else if (covs == "true_Si") {  
    # Probit model includes C and S: C + S + C:S
    svy_data <- data.frame(s = factor(s_all),
                           x = x_mat,
                           y = y_all, 
                           wts = w_all,
                           c = factor(analysis$c_all),
                           clus = clus_id_all)
    svydes <- svydesign(ids = ~clus, strata = ~s, weights = ~wts, data = svy_data)
    fit <- svyglm(y ~ c * s, design = svydes, family = quasibinomial(link = "probit"))                     
  } else if (covs == "additional") {
    # Probit model includes C, S, A (binary), and B (continuous)
    # C + S + A + B + C:S + C:A + C:B
    svy_data <- data.frame(s = factor(s_all),
                           x = x_mat,
                           y = y_all, 
                           wts = w_all,
                           c = factor(analysis$c_all),
                           a = factor(a_all),
                           b = b_all,
                           clus = clus_id_all)
    svydes <- svydesign(ids = ~clus, strata = ~s, weights = ~wts, data = svy_data)
    fit <- svyglm(y ~ c + s + a + b + c:s + c:a + c:b, 
                  design = svydes, family = quasibinomial(link = "probit")) 
  } 
  coefs <- fit$coefficients
  ci <- confint(fit)
  
  # Convert format to match WSOLCA and SOLCA
  xi_med <- xi_med_lb <- xi_med_ub <- matrix(NA, nrow = analysis$K_red, ncol = q)
  if (is.null(covs)) {
    # Probit model only includes latent class C
    for (k in 1:analysis$K_red) {
      xi_med[k, 1] <- coefs[1] + (k != 1) * coefs[q + (k-1)] 
      xi_med_lb[k, 1] <- ci[1, 1] + (k != 1) * ci[q + (k-1), 1] 
      xi_med_ub[k, 1] <- ci[1, 2] + (k != 1) * ci[q + (k-1), 2]
    }
  } else {
    # Probit model includes C, S, and potentially also A and B
    xi_med[1, ] <- coefs[c(1, (analysis$K_red + 1:(q-1)))]
    xi_med_lb[1, ] <- ci[c(1, (analysis$K_red + 1:(q-1))), 1]
    xi_med_ub[1, ] <- ci[c(1, (analysis$K_red + 1:(q-1))), 2]
    for (k in 2:analysis$K_red) {
      xi_med[k, ] <- coefs[c(k, 
        (k + (q-1)) + (analysis$K_red-1) * (1:(q-1)))] + xi_med[1, ]
      xi_med_lb[k, ] <- ci[c(k, 
        (k + (q-1)) + (analysis$K_red-1) * (1:(q-1))), 1] + xi_med_lb[1, ]
      xi_med_ub[k, ] <- ci[c(k, 
        (k + (q-1)) + (analysis$K_red-1) * (1:(q-1))), 2] + xi_med_ub[1, ]
    }
  }
    
  analysis$xi_med <- xi_med
  analysis$xi_med_lb <- xi_med_lb
  analysis$xi_med_ub <- xi_med_ub
  
  runtime <- Sys.time() - start_time
  
  #================= Save and return output ====================================
  res <- list(analysis = analysis, runtime = runtime, 
              data_vars = data_vars, V = V, MCMC_out = MCMC_out, 
              post_MCMC_out = post_MCMC_out, K_MCMC = adapt_MCMC$K_MCMC)
  if (save_res) {
    save(res, file = res_path)
  }
  return(res)
}



#===================== RUN MAIN WOLCA FUNCTION =================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/SWOLCA/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/SWOLCA/"
data_dir <- "Data/July6/"
res_dir <- "Results/July6/"
model_dir <- "Model_Code/"
model <- "wOFMM"

# # Testing code
# scen_samp <- 111211
# iter_pop <- 1
# samp_n <- 1
# n_runs <- 100
# burn <- 50
# thin <- 5
# save_res <- FALSE
# covs <- "true_Si"

# Define paths
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".RData")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_scen", scen_samp, 
                     "_samp", samp_n, ".RData")  # Output file
res_path <- paste0(wd, res_dir, model, "_results_wt_scen", scen_samp, 
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
  print(paste0("Running WOLCA_main for scenario ", scen_samp, ' iter ', 
               iter_pop,' samp ', samp_n))
  results <- WOLCA_main_Rcpp(data_path = data_path, adapt_path = adapt_path,
                             res_path = res_path,
                             save_res = save_res, n_runs = n_runs, 
                             burn = burn, thin = thin, covs = covs)
  print(paste0("Runtime: ", results$runtime))
}

