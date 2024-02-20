#' @name baysc-package
#' @aliases baysc-package baysc
# #' @docType package
#' @title BAYesian Survey Clustering
#' 
#' @description
#' `baysc` is an R package for running Bayesian clustering methods on survey data.
#' 
#' @details
#' Use the `swolca()` function to run a supervised weighted overfitted latent 
#' class analysis (SWOLCA), proposed in Wu et al. (2023). 
#' 
#' Use the `solca()` function to run an unweighted supervised overfitted latent 
#' class analysis (SOLCA) that does not account for survey design.
#' 
#' Use the `wolca()` function to run a two-step analysis where the first step is
#' an unsupervised weighted overfitted latent class analysis (WOLCA) and the 
#' second step is a frequentist survey-weighted probit regression model using 
#' the `survey` package.
#' 
#' A simulated dataset named `sim_data` is also provided to help with 
#' familiarization of the package.
#' 
#' @references Wu, S. M., Williams, M. R., Savitsky, T. D., & Stephenson, B. J. 
#' (2023). Derivation of outcome-dependent dietary patterns for low-income women 
#' obtained from survey data using a Supervised Weighted Overfitted Latent Class 
#' Analysis. arXiv preprint arXiv:2310.01575.
#' 
#' @keywords package
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' 
#' @examples
#' \dontrun{
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
#' stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
#' sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
#' n <- dim(x_mat)[1]        # Number of individuals
#' 
#' # Probit model only includes latent class
#' V <- matrix(1, nrow = n) # Regression design matrix without class assignment
#'
#' # Stan model
#' mod_stan <- stanmodels$SWOLCA_main
#' 
#' # Run swolca
#' res_swolca <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
#' cluster_id = cluster_id, stratum_id = stratum_id, V = V, adapt_seed = 1,
#' n_runs = 50, burn = 25, thin = 1, mod_stan = mod_stan, save_res = FALSE)
#' 
#' # Run solca
#' res_solca <- solca(x_mat = x_mat, y_all = y_all, V = V, adapt_seed = 1, 
#' n_runs = 50, burn = 25, thin = 1, save_res = FALSE)
#' 
#' # Run wolca
#' res_wolca <- wolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
#' cluster_id = cluster_id, stratum_id = stratum_id, V = V, glm_form = glm_form, 
#' adapt_seed = 1, n_runs = 50, burn = 25, thin = 1, mod_stan = mod_stan, 
#' save_res = FALSE)
#' }
#' @keywords internal
"_PACKAGE"