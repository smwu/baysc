#' Fit survey-weighted probit model for WOLCA
#'
#' @description
#' `wolca_svyglm` relates latent classes patterns derived from [wolca()] to a 
#' binary outcome by fitting a survey-weighted probit model.
#'
#' @inheritParams swolca
#' @param res An object of class `"wolca"`, resulting from a call to [wolca()] 
#' or [swolca_var_adjust()]
#' @param ci_level Confidence interval level for probit regression coefficient 
#' estimates. Default is `0.95`.
#' 
#' @details
#' `wolca_svyglm` is the second step of a two-step approach that runs an 
#' unsupervised WOLCA in the first step to derive latent class patterns and 
#' subsequently treats the class assignments as fixed and includes them as 
#' covariates in a frequentist survey-weighted probit regression model that uses 
#' an asymptotic sandwich estimator for variance estimation. 
#' 
#' `wolca_svyglm` specifies survey design and fits a survey-weighted probit 
#' regression model according to the formula specified in `glm_form` using the 
#' `svyglm()` function from the `survey` package (Lumley, 2023). Regression 
#' coefficients and their confidence intervals are obtained from the `svyglm()` 
#' output. If the residual degrees of freedom is negative, a Wald confidence 
#' interval is manually calculated using a t-distribution with degrees of 
#' freedom from the survey design, and the `fit_summary` output also uses the 
#' degrees of freedom from the survey design to obtain p-values. This can happen 
#' when there are few clusters per stratum and a large number of covariates. 
#' 
#' The point and interval estimates are then converted into mixture 
#' reference coding format to match the output format from [swolca()]. `V_data` 
#' includes all covariates to include in the probit regression other than latent 
#' class. 
#' 
#' To save results, set `save_res = TRUE` (default) and `save_path` to a string
#' that specifies both the location and the beginning of the file name 
#' (e.g., "~/Documents/run"). The file name will have "_wolca_results.RData" 
#' appended to it, overwriting the unsupervised results if the file names are the 
#' same.
#' 
#' @return
#' Returns an object `res` of class `"wolca"`, which includes all outputs from 
#' [wolca()] or [wolca_var_adjust()] as well as additional list `estimates_syglm` 
#' containing the following objects:
#' \describe{
#'   \item{\code{xi_est}}{Matrix of estimates for xi. (K_red)xQ}
#'   \item{\code{xi_est_lb}}{Matrix of confidence interval lower bound estimates for xi. (K_red)xQ}
#'   \item{\code{xi_est_ub}}{Matrix of confidence interval upper bound estimates for xi. (K_red)xQ}
#'   \item{\code{fit}}{`svyglm` class object with output from the `svyglm` regression model}
#'   \item{\code{fit_summary}}{`summary.svyglm` class object with output from 
#' the `summary()` function. If the residual degrees of freedom is negative, 
#' the degrees of freedom from the survey design is used.}
#' }
#' List `data_vars` of `res` is updated to contain the following additional objects:
#' \describe{
#'   \item{\code{Q}}{Number of regression covariates excluding class assignment.}
#'   \item{\code{y_all}}{Vector of binary outcomes; nx1.}
#'   \item{\code{V_data}}{Dataframe of additional regression covariates; nxQ or `NULL`.}
#'   \item{\code{V}}{Regression design matrix without class assignment; nxQ.}
#'   \item{\code{glm_form}}{String specifying formula for probit regression, 
#'   excluding outcome and latent class.}
#'   \item{\code{ci_level}}{Confidence interval level for probit regression 
#'   coefficient estimates.}
#' }
#' The `runtime` output for `res` is also updated to include the runtime for the 
#' probit regression in addition to the runtime for the main `wolca()` model.
#' 
#' If `save_res = TRUE` (default), the updated `res` object is saved as 
#' `[save_path]_wolca_results.RData`, overwriting the unsupervised results if 
#' the file names are the same. 
#'
#' @seealso [wolca()] 
#' @importFrom survey svydesign svyglm degf
#' @importFrom stats confint as.formula quasibinomial terms as.formula
#' @export
#' 
#' @references Lumley T (2023). “survey: analysis of complex survey samples.” R 
#' package version 4.2.
#'
#' @examples
#' \dontrun{   
#' # Load data and obtain relevant variables
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
#' # Survey-weighted regression formula
#' glm_form <- "~ 1"
#' 
#' # Run wolca
#' res <- wolca(x_mat = x_mat, sampling_wt = sampling_wt, 
#'              cluster_id = cluster_id, stratum_id = stratum_id, 
#'              run_sampler = "both", adapt_seed = 1, n_runs = 50, burn = 25, 
#'              thin = 1, save_res = FALSE)
#'        
#' # Apply variance adjustment to posterior estimates
#' res_adjust <- wolca_var_adjust(res = res, num_reps = 100, save_res = FALSE, 
#'                                adjust_seed = 1) 
#' 
#' # Run weighted outcome regression model
#' res_svyglm <- wolca_svyglm(res = res_adjust, y_all = y_all, 
#'                            glm_form = glm_form, ci_level = 0.95, 
#'                            V_data = V_data, save_res = FALSE)
#' }
wolca_svyglm <- function(res, y_all, V_data = NULL, glm_form, ci_level = 0.95, 
                         save_res = TRUE, save_path = NULL) {
  
  #============== Catch errors and initialize variables ========================
  # Check object class and estimates
  if (!inherits(res, "wolca")) {
    stop("res must be an object of class `wolca`, resulting from a call to the 
         `wolca()` function that includes results from the fixed sampler")
  } else if (is.null(res$estimates)) {
    stop("res must include results from the fixed sampler in the `wolca()` function")
  }
  
  # Begin runtime tracker
  start_time <- Sys.time()
  
  # Extract data elements into the global environment
  stratum_id <- res$data_vars$stratum_id
  cluster_id <- res$data_vars$cluster_id
  x_mat <- res$data_vars$x_mat
  w_all <- res$data_vars$w_all
  
  # If no additional covariates, set V_data to be a column of all ones
  if (is.null(V_data)) {
    V_data <- as.data.frame(matrix(1, nrow = length(w_all)))
  }
  
  # Catch errors
  catch_errors(x_mat = res$data_vars$x_mat, y_all = y_all, V_data = V_data, 
               glm_form = glm_form, save_res = save_res, save_path = save_path)
  
  # Obtain probit regression design matrix without class assignment
  V <- model.matrix(as.formula(glm_form), data = V_data)
  # Number of regression covariates excluding class assignment 
  Q <- ncol(V)  
  
  # Set pointer to adjusted or unadjusted estimates
  if (!is.null(res$estimates_adjust)) {
    # Adjusted estimates
    estimates <- res$estimates_adjust
  } else {
    # Unadjusted estimates
    estimates <- res$estimates
  }
  
  #============== Create survey design =========================================
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
  
  #============== Fit probit model for the outcome =============================
  # Add outcome and latent class main and interaction terms to formula
  terms <- labels(stats::terms(stats::as.formula(glm_form)))
  if (length(terms) > 0) {
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
  fit_summary <- summary(fit)
  # Obtain coefficients and confidence interval
  coefs <- fit$coefficients
  ci <- stats::confint(fit, level = ci_level)
  
  # Residual df is the design df (the number of PSUs minus the number of strata) 
  # minus the number of predictors, which can easily become negative when you 
  # have only a few large clusters per stratum.
  if (all(is.na(ci))) {
    # If zero/negative residual df, manually calculate the Wald confidence interval 
    # using a t-distribution with degrees of freedom from the survey design. 
    # Can also use survey degf for the `summary()` function and get p-values
    # If no covariates at the cluster level, there is a reasonable argument that 
    # the regression doesn't use up degrees of freedom, so can use svy design df 
    ci <- manual_CI(model_object = fit, svy_df = survey::degf(svydes), 
                    ci = ci_level)[, -1]
    fit_summary <- summary(fit, df = survey::degf(svydes))
    warning("Residual degrees of freedom negative. Degrees of freedom from the survey design used instead for calculating CI and p-values.")
  }
  
  # Convert format to mixture reference to match SWOLCA and SOLCA
  xi_list <- convert_ref_to_mix(K = length(estimates$pi_med), Q = Q, 
                                est_beta = coefs, ci_beta = ci)
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  # Add probit regression runtime to overall runtime
  sum_runtime <- runtime + res$runtime
  res$runtime <- sum_runtime
  
  # Add probit regression estimates to output
  res$estimates_svyglm <- list()
  res$estimates_svyglm$xi_est <- xi_list$est_xi
  res$estimates_svyglm$xi_est_lb <- xi_list$est_xi_lb
  res$estimates_svyglm$xi_est_ub <- xi_list$est_xi_ub
  res$estimates_svyglm$fit <- fit
  res$estimates_svyglm$fit_summary <- fit_summary
  
  # Add outcome data to output
  res$data_vars$Q = Q
  res$data_vars$y_all = y_all
  res$data_vars$V_data = V_data
  res$data_vars$V = V
  res$data_vars$glm_form = glm_form
  res$data_vars$ci_level = ci_level
  
  class(res) <- "wolca"
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_wolca_results.RData"))
  }
  
  return(res)
}
