#' Fit probit model for WOLCA
#'
#' @description
#' `fit_probit_wolca` uses `svyglm` to fit a survey-weight probit model in a 
#' two-step model where the first step derived latent classes using an 
#' unsupervised WOLCA 
#'
#' @inheritParams wolca
#' @param estimates Output from `get_estimates_wolca()` containing `K_red`, 
#' `pi_red`, `theta_red`, `pi_med`, `theta_med`, `c_all`, `pred_class_probs`
#' @param w_all Weights normalized to sum to n. nx1
#' @param q Number of regression covariates excluding class assignment
#' 
#' @details
#' Specifies survey design and fits a survey-weighted probit regression model
#' according to the formula specified in `glm_form`. Regression coefficients and 
#' their confidence intervals are obtained from the `svyglm()` output. If the 
#' residual degrees of freedom is less than 1, a Wald confidence interval is 
#' manually calculated using a t-distribution with degrees of freedom from the 
#' survey design. The point and interval estimates are then converted into the 
#' factor reference coding format to match the output from `swolca()` and `solca()`. 
#' 
#' @return
#' Returns updated list `estimates` containing the following additional objects:
#' \describe{
#'   \item{\code{xi_est}}{Matrix of estimates for xi. (K_red)xq}
#'   \item{\code{xi_est_lb}}{Matrix of confidence interval lower bound estimates for xi. (K_red)xq}
#'   \item{\code{xi_est_ub}}{Matrix of confidence interval upper bound estimates for xi. (K_red)xq}
#'   \item{\code{fit}}{`svyglm` class object with output from the `svyglm` regression model}
#' }
#'
#' @seealso [run_MCMC_Rcpp_wolca()] [post_process_wolca()] 
#' [get_estimates_wolca()] [wolca()] 
#' @importFrom survey svydesign svyglm degf
#' @importFrom stats confint as.formula quasibinomial terms as.formula
#' @export
#'
#' @examples
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
#' burn = 25, thin = 5, K = K, J = J, R = R, n = n, w_all = w_all, x_mat = x_mat, 
#' alpha = alpha, eta = eta)
#' 
#' # Then run post-process relabeling
#' post_MCMC_out <- post_process_wolca(MCMC_out = MCMC_out, J = J, R = R,
#' class_cutoff = 0.05)
#'
#' # Then obtain posterior estimates for WOLCA
#' estimates <- get_estimates_wolca(MCMC_out = MCMC_out, 
#' post_MCMC_out = post_MCMC_out, n = n, J = J, x_mat = x_mat)
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
#' stratum_id = stratum_id, cluster_id = cluster_id, x_mat = x_mat, 
#' y_all = y_all, w_all = w_all, V_data = V_data, q = q)
#' 
fit_probit_wolca <- function(estimates, glm_form, stratum_id, cluster_id, 
                             x_mat, y_all, w_all, ci_level = 0.95, V_data, q) {
  
  # Create survey design
  if (!is.null(stratum_id)) {  # Include stratifying variable
    # Survey data frame for specifying survey design
    svy_data <- data.frame(stratum_id = factor(stratum_id), 
                           cluster_id = factor(cluster_id),
                           x_mat = x_mat, y_all = y_all, w_all = w_all)
    # Add latent class assignment variable to survey data
    svy_data$c_all <- factor(estimates$c_all)
    # Add additional covariates
    svy_data <- cbind(svy_data, V_data)
    # Specify survey design
    svydes <- survey::svydesign(ids = ~cluster_id, strata = ~stratum_id, 
                                weights = ~w_all, data = svy_data)
  } else { # No stratifying variable
    # Survey data frame for specifying survey design
    svy_data <- data.frame(cluster_id = cluster_id, x_mat = x_mat, 
                           y_all = y_all, w_all = w_all)
    # Add latent class assignment variable to survey data
    svy_data$c_all <- factor(estimates$c_all)
    # Add additional covariates
    svy_data <- cbind(svy_data, V_data)
    # Specify survey design
    svydes <- survey::svydesign(ids = ~cluster_id, weights = ~w_all, 
                                data = svy_data)    
  }
  
  # Add outcome and latent class main and interaction terms to formula
  terms <- labels(stats::terms(stats::as.formula(glm_form)))
  if (length(terms) > 1) {
    full_glm_form <- paste0("y_all ~ ", 
                            paste0("c_all * ", terms, collapse = " + ")) 
  } else {
    full_glm_form <- paste0("y_all ~ c_all") 
  }
  
  # If only one latent class, cannot have latent class as a covariate
  if (length(levels(svy_data$c_all)) == 1) {
    stop("Only one latent class found. Cannot use latent class as a covariate")
  }
  
  # Fit probit model according to specified formula
  fit <- survey::svyglm(formula = stats::as.formula(full_glm_form), 
                        design = svydes, 
                        family = stats::quasibinomial(link = "probit"))
  # Obtain coefficients and confidence interval
  coefs <- fit$coefficients
  ci <- stats::confint(fit, level = ci_level)
  # If zero/negative residual df, manually calculate the Wald confidence interval 
  # using a t-distribution with degrees of freedom from the survey design. 
  # Best if no cluster-level covariates in the regression model
  if (all(is.na(ci))) {
    ci <- manual_CI(model_object = fit, svy_df = survey::degf(svydes), 
                    ci = ci_level)[, -1]
  }
  
  # Convert format to match SWOLCA and SOLCA
  xi_list <- convert_ref_to_mix(K = estimates$K_red, q = q, est_beta = coefs,
                               ci_beta = ci)

  # Return output with probit regression estimates
  estimates$xi_est <- xi_list$est_xi
  estimates$xi_est_lb <- xi_list$est_xi_lb
  estimates$xi_est_ub <- xi_list$est_xi_ub
  estimates$fit <- fit
  
  return(estimates)
}
