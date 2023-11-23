#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' cluster_id <- data_vars$cluster_id  # Cluster indicators, nx1
#' stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
#' sampling_wt <- data_vars$sample_wt
#' 
#' # Obtain dimensions
#' n <- dim(x_mat)[1]        # Number of individuals
#' J <- dim(x_mat)[2]        # Number of exposure items
#' R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
#'              function(x) length(unique(x)))
#' R <- max(R_j)             # Maximum number of exposure categories across items
#' # Obtain normalized weights
#' kappa <- sum(sampling_wt) / n   # Weights norm. constant
#' w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
#' 
#' # Set hyperparameters for fixed sampler
#' K <- 3
#' alpha <- rep(1, K) / K
#' eta <- matrix(0.01, nrow = J, ncol = R)
#' for (j in 1:J) {
#'   eta[j, 1:R_j[j]] <- rep(1, R_j[j])
#' }
#' 
#' # First initialize OLCA params
#' OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)
#' 
#' # Then run MCMC sampling
#' MCMC_out <- run_MCMC_Rcpp_wolca(OLCA_params = OLCA_params, n_runs = 50,
#'                                 burn = 25, thin = 5, K = K, J = J, R = R, n = n, w_all = w_all, x_mat = x_mat,
#'                                 alpha = alpha, eta = eta)
#' 
#' # Then run post-process relabeling
#' post_MCMC_out <- post_process_wolca(MCMC_out = MCMC_out, J = J, R = R,
#'                                     class_cutoff = 0.05)
#' 
#' # Then obtain posterior estimates for WOLCA
#' estimates <- get_estimates_wolca(MCMC_out = MCMC_out,
#'                                  post_MCMC_out = post_MCMC_out, n = n, J = J, x_mat = x_mat)
#' 
#' # Probit model only includes latent class
#' V_data <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
#' # Survey-weighted regression formula
#' glm_form <- "~ 1"
#' # Obtain probit regression design matrix without class assignment
#' V <- model.matrix(as.formula(glm_form), data = V_data)
#' # Number of regression covariates excluding class assignment
#' q <- ncol(V)
#' 
#' # Finally run weighted probit regression model
#' estimates <- fit_probit_wolca(estimates = estimates, glm_form = glm_form,
#'                               stratum_id = stratum_id, cluster_id = cluster_id, x_mat = x_mat,
#'                               y_all = y_all, w_all = w_all, V_data = V_data, q = q)
#' 
#' 
#' # Probit model only includes latent class
#' V_data <- data.frame(stratum_id = stratum_id) # Additional regression covariates
#' # Survey-weighted regression formula
#' glm_form <- "~ stratum_id"
#' 
#' # Run wolca
#' res <- wolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
#'              cluster_id = cluster_id, stratum_id = stratum_id, V_data = V_data,
#'              run_sampler = "both", glm_form = glm_form, adapt_seed = 1,
#'              n_runs = 50, burn = 25, thin = 1, save_res = FALSE)
#' 
#' res2 <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
#'                cluster_id = cluster_id, stratum_id = stratum_id, V_data = V_data,
#'                run_sampler = "both", glm_form = glm_form, adapt_seed = 1,
#'                n_runs = 50, burn = 25, thin = 1, save_res = FALSE)
#' get_regr_coefs(res_wolca)
#' get_regr_coefs(res_swolca, digits = 3)
#' get_regr_coefs(res2)
#' 
#' #' @param format String indicating format for the kable output. Possible values
#' #' are `"latex"` (default), `"html"`, `"pipe"`, `"simple"`, `"rst"`, `"jira"`,
#' #' and `"org"`. See documentation for `kable()` for more information.