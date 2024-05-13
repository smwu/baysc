#===================================================
## Helper functions for SWOLCA and WOLCA
## Programmer: SM Wu   
## Data: Simulations and application   
#===================================================


#' Get mode
#' 
#' `get_mode` is a helper function that obtains the modal value given an input 
#' vector.
#' @param v Input vector
#' @return Outputs most common value found in input vector `v`
#' @keywords internal
#' @export
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#' Get manual confidence interval
#' 
#' `manual_CI` manually calculates a Wald confidence interval using a t-dist
#' with df from the survey design. Used for WOLCA in the situation where svyglm 
#' produces negative residual df, calculated as design df plus one, minus the 
#' number of parameters estimated. Best if no cluster-level covariates in the 
#' regression model
#' @param model_object svyglm model fit object
#' @param svy_df survey design df
#' @param  ci confidence interval level
#' @return Outputs dataframe of confidence interval for all coefficients
#' @importFrom stats coef qt
#' @keywords internal
#' @export
manual_CI <- function(model_object, svy_df, ci = 0.95){
  a <- stats::coef(summary(model_object))
  mult <- stats::qt((1 + ci) / 2, df = svy_df)
  restab <- with(as.data.frame(a),
                 cbind(est = Estimate,
                       lwr =  Estimate - mult*`Std. Error`,
                       upr = Estimate + mult*`Std. Error`))
  rownames(restab) <- rownames(a)
  return(data.frame(restab))
}


#' Catch errors in input variables for external functions
#' 
#' @description
#' Catch input errors in variables necessary for package functions. All 
#' parameters are set to `NULL` by default so that error checks 
#' are only performed on relevant variables.
#' 
#' @inheritParams swolca
#' @param model String specifying which model is used. Must be one of `swolca` 
#' (default) or `wolca`
#' @return The function stops and an error message is displayed if the input 
#' variables are not acceptable
#' @details All parameters are set to `NULL` by default so that error checks 
#' are only performed on relevant variables.
#' 
#' @importFrom stats terms as.formula
#' @importFrom stringr str_detect
#' @keywords internal
#' @export
catch_errors <- function(x_mat = NULL, y_all = NULL, sampling_wt = NULL, 
                         cluster_id = NULL, stratum_id = NULL, V_data = NULL,
                         run_sampler = NULL, glm_form = NULL, Q = NULL,
                         K_max = NULL, adapt_seed = NULL, fixed_seed = NULL,
                         class_cutoff = NULL, alpha_adapt = NULL, 
                         eta_adapt = NULL, mu0_adapt = NULL, 
                         Sig0_adapt = NULL, K_fixed = NULL, alpha_fixed = NULL, 
                         eta_fixed = NULL, mu0_fixed = NULL, Sig0_fixed = NULL,
                         n_runs = NULL, burn = NULL, thin = NULL, update = NULL,
                         adjust_var = NULL, num_reps = NULL,
                         save_res = NULL, save_path = NULL) {
  if (is.null(x_mat)) {
    stop("need to specify exposure matrix")
  } else {
    
    # Check type for x_mat
    if (!is.matrix(x_mat)) {
      stop("x_mat must be a matrix")
    }
    
    # Obtain dimensions
    n <- dim(x_mat)[1]        # Number of individuals
    J <- dim(x_mat)[2]        # Number of exposure items
    R <- max(apply(x_mat, 2,  # Number of exposure categories
                   function(x) length(unique(x))))  
    
    # Check sampler specification
    if (!is.null(run_sampler)) {
      if (!(run_sampler %in% c("both", "adapt", "fixed"))) {
        stop("run_sampler must be one of `both`, `adapt`, or `fixed`")
      }
      if (run_sampler == "fixed") {
        if (is.null(K_fixed)) {
          stop("K_fixed must be specified")
        }
      }
    }
    
    # Check regression formula
    if (!is.null(glm_form)) {
      if (substring(glm_form, 1, 1) != "~") {
        stop("glm_form must be a string starting with '~' that specifies a valid 
             formula")
      }
      regr_vars <- labels(stats::terms(stats::as.formula(glm_form)))   
      if (any(stringr::str_detect(regr_vars, "c_all"))) {
        stop("please exclude latent class assignment, c_all, from glm_form, 
              as it is already assumed to be included")
      }
      # Extract additional covariates, not including interaction terms
      regr_vars <- regr_vars[!(stringr::str_detect(regr_vars, ":"))]
      if (!(all(regr_vars %in% colnames(V_data)))) {
        stop("all variables in glm_form must exist in V_data")
      }
    }
    
    # Check same number of individuals for x and y
    if (!is.null(y_all)) {
      if (n != length(y_all)) {
        stop("number of rows in x_mat must match length of y_all")
      }
      if (!is.vector(y_all)) {
        stop("y_all must be a vector")
      }
    }
    
    # Check same number of individuals for x and sampling weights
    if (!is.null(sampling_wt)) {
      if (!is.vector(sampling_wt)) {
        stop("sampling_wt must be a vector")
      }
      if (n != length(sampling_wt)) {
        stop("number of rows in x_mat must match length of sampling_wt")
      }
    }
    
    # If no clustering, assign each individual to their own cluster. Else, check
    # same number of individuals for x and clusters
    if (is.null(cluster_id)) {
      cluster_id <- 1:n
    } else {
      if (!is.vector(cluster_id)) {
        stop("cluster_id must be a vector")
      }
      if (n != length(cluster_id)) {
        stop("number of rows in x_mat must match length of cluster_id")
      }
    }
    
    # Check same number of individuals for x and strata
    if (!is.null(stratum_id)) {
      if (!is.vector(stratum_id)) {
        stop("stratum_id must be a vector")
      }
      if (n != length(stratum_id)) {
        stop("number of rows in x_mat must match length of stratum_id")
      }
    }
    
    # Check same number of individuals for x and V_data
    if (!is.null(V_data)) {
      if (n != nrow(V_data)) {
        stop("number of rows in x_mat must match number of rows in V_data")
      }
      if (class(V_data)[1] != "data.frame") {
        stop("V_data must be a dataframe")
      }
    }
    
    # Check cutoff is between 0 and 1
    if (!is.null(class_cutoff)) {
      if (class_cutoff <= 0 | class_cutoff >= 1) {
        stop("class_cutoff must be a proportion in (0,1)")
      }
    }
    
    # Check seeds
    if (!is.null(adapt_seed)) {
      if (!is.numeric(adapt_seed)) {
        stop("adapt_seed must be numeric")
      }
    }
    if (!is.null(fixed_seed)) {
      if (!is.numeric(fixed_seed)) {
        stop("fixed_seed must be numeric")
      }
    }
    
    # Check hyperparameter dimensions for adaptive sampler
    if (!is.null(alpha_adapt)) {
      if (length(alpha_adapt) != K_max) {
        stop("length of alpha_adapt must be the same as K_max")
      }
    }
    if (!is.null(eta_adapt)) {
      if ((nrow(eta_adapt) != J) | (ncol(eta_adapt) != R)) {
        stop("eta_adapt must be a matrix with J rows and R columns, where J is
             the number of exposure items and R is the maximum number of 
             exposure categories")
      }
    }  
    if (!is.null(mu0_adapt)) {
      if (length(mu0_adapt) != K_max | !is.list(mu0_adapt) | 
          !(all(lapply(mu0_adapt, length) == Q))) {
        stop("mu0_adapt must be a list of length K_max where each element is a 
           vector of length Q (number of regression covariates excluding latent class)")
      }
    }
    if (!is.null(Sig0_adapt)) {
      if (length(Sig0_adapt) != K_max | !is.list(Sig0_adapt) | 
          !(all(lapply(Sig0_adapt, nrow) == Q)) | 
          !(all(lapply(Sig0_adapt, ncol) == Q))) {
        stop("Sig0_adapt must be a list of length K_max where each element is a 
            QxQ matrix, where Q is the number of regression covariates excluding 
           latent class)")
      }
    }
    
    # Check number of classes
    if (!is.null(K_max)) {
      if (K_max < 1) {
        stop("Maximum number of classes must be at least 1")
      }
    }
    if (!is.null(K_fixed)) {
      if (K_fixed < 1) {
        stop("Maximum number of classes must be at least 1")
      }
      # Check hyperparameter dimensions for fixed sampler
      if (!is.null(alpha_fixed)) {
        if (length(alpha_fixed) != K_fixed) {
          stop("length of alpha_fixed must be the same as K_fixed")
        }
      }
      if (!is.null(eta_fixed)) {
        if ((nrow(eta_fixed) != J) | (ncol(eta_fixed) != R)) {
          stop("eta_fixed must be a matrix with J rows and R columns, where J is
             the number of exposure items and R is the maximum number of 
             exposure categories")
        }
        if (any(eta_fixed == 0)) {
          warning("eta_fixed has 0 values and may result in rank-difficiency issues
                during the Hessian calculation in the var_adjust() function")
        }
      }
      if (!is.null(mu0_fixed)) {
        if (length(mu0_fixed) != K_fixed | !is.list(mu0_fixed) | 
            !(all(lapply(mu0_fixed, length) == Q))) {
          stop("mu0_fixed must be a list of length K_fixed where each element is a 
           vector of length Q (number of regression covariates excluding latent class)")
        }
      }
      if (!is.null(Sig0_fixed)) {
        if (length(Sig0_fixed) != K_fixed | !is.list(Sig0_fixed) | 
            !(all(lapply(Sig0_fixed, nrow) == Q)) | 
            !(all(lapply(Sig0_fixed, ncol) == Q))) {
          stop("Sig0_fixed must be a list of length K_fixed where each element is a 
            QxQ matrix, where Q is the number of regression covariates excluding 
           latent class)")
        }
      }
    } else {
      if (any(!is.null(c(alpha_fixed, eta_fixed, mu0_fixed, Sig0_fixed)))) {
        stop("K_fixed must be specified alongside the priors for the fixed sampler")
      }
    }
    
    # Check MCMC parameters
    if (!all(is.null(c(n_runs, burn, thin, update)))) {
      if (!all(c(n_runs, burn, thin, update) %% 1 == 0) | 
          !all(c(n_runs, burn, thin, update) >= 0)) {
        stop("n_runs, burn, thin, and update must be whole numbers")
      }
      if (burn > n_runs) {
        stop("n_runs must be larger than burn")
      }
    }
    
    # Check variance adjustment parameters
    if (!is.null(adjust_var)) {
      if (!is.logical(adjust_var)) {
        stop("adjust_var must be a boolean specifying if the post-processing 
        variance adjustment for accurate interval coverage should be applied")
      } else if (adjust_var) {
        if (!is.integer(num_reps) | num_reps < 1) {
          stop("num_reps must be a positive integer, recommended to be at least
          50. More replicates will lead to more accurate results but will take 
          longer to run.")
        }
      }
    } 
    
    # Check saving parameters
    if (!is.null(save_res)) {
      if (!is.logical(save_res)) {
        stop("save_res must be a boolean specifying if results should be saved")
      }
      if (save_res) {
        if (is.null(save_path) | !is.character(save_path)) {
          stop("save_path must be a string specifying a path and file name, such as '~/Documents/run'")
        } else {
          last_slash_ind <- regexpr("\\/[^\\/]*$", save_path)
          if (!dir.exists(substr(save_path, start = 1, stop = last_slash_ind))) {
            stop("directory specified in save_path does not exist")
          }
          if (last_slash_ind == length(save_path)) {
            stop("please append the start of a file name to the end of save_path. 
            For example, '~/Documents/run' can produce a saved file named 
            'run_swolca_results.RData'")
          }
        }
      }
    }
  }
}


#' Convert from reference cell coding to mixture reference coding 
#'
#' `convert_ref_to_mix` converts from reference cell coding to a combination of 
#' factor variable and reference cell coding, referred to as mixture reference 
#' coding.
#' 
#' @param K Number of latent classes
#' @param Q Number of regression covariates excluding class assignment
#' @param est_beta Vector of probit regression coefficients in reference cell 
#' coding. (K*Q)x1. Order of betas must correspond to a formula with `c_all` as
#' the first covariate and all interaction terms involving `c_all` present.
#' @param ci_beta Matrix of interval estimates for probit regression 
#' coefficients in reference cell coding. (K*Q)x2, where the first column is 
#' the lower bound and the second column is the upper bound. Set to `NULL` 
#' (default), if no interval estimate conversions are necessary.
#' 
#' @return Outputs list `xi_list` containing `est_xi`, a KxQ matrix of 
#' regression estimates in mixture reference coding. If `ci_beta` is not `NULL`,
#' `xi_list` also contains `est_xi_lb` and `est_xi_ub`, which are KxQ matrices
#' of the lower and upper bound interval estimates, respectively, in mixture 
#' reference coding.
#' 
#' @keywords internal
#' @seealso [wolca_svyglm()]
#' @export
convert_ref_to_mix <- function(K, Q, est_beta, ci_beta = NULL) {
  est_xi <- matrix(NA, nrow = K, ncol = Q)
  if (!is.null(ci_beta)) {
    est_xi_lb <- est_xi_ub <- matrix(NA, nrow = K, ncol = Q)
  }
  
  # If only latent class as a covariate (no interactions)
  if (Q == 1) {
    est_int_k1 <- NULL   # baseline class
    est_int_koth <- NULL # additional classes
  # Position of interaction terms for additional covariates
  } else if (Q > 1) {
    est_int_k1 <- K + 1:(Q-1)       # baseline class
    est_int_koth <- (K - 1) * (1:(Q-1)) # additional classes
  }
  
  # Get estimates for first class level, including all interactions
  est_xi[1, ] <- est_beta[c(1, est_int_k1)]
  if (!is.null(ci_beta)) {
    est_xi_lb[1, ] <- ci_beta[c(1, est_int_k1), 1]
    est_xi_ub[1, ] <- ci_beta[c(1, est_int_k1), 2]
  }
  # Get estimates for additional class levels, including all interactions
  for (k in 2:K) {
    est_xi[k, ] <- est_beta[c(k, (k + (Q-1)) + est_int_koth)] + est_xi[1, ]
    if (!is.null(ci_beta)) {
      est_xi_lb[k, ] <- ci_beta[c(k, (k + (Q-1)) + est_int_koth), 1] + est_xi_lb[1, ]
      est_xi_ub[k, ] <- ci_beta[c(k, (k + (Q-1)) + est_int_koth), 2] + est_xi_ub[1, ]
    }
  }
  
  # Return estimates and credible intervals in mixture reference coding
  xi_list <- list(est_xi = est_xi)
  if (!is.null(ci_beta)) {
    xi_list$est_xi_lb = est_xi_lb
    xi_list$est_xi_ub = est_xi_ub
  }
  return(xi_list)
      # beta_comb <- beta_ref
      # for (i in 2:nrow(beta_ref)) {
      #   beta_comb[i, ] <- beta_ref[i, ] - beta_ref[1, ]
      # }
      # return(beta_comb)
}


#' Convert from mixture reference coding to reference cell coding
#' 
#' @description
#' Convert regression estimates \eqn{\xi} from mixture reference coding, a 
#' combination of factor variable and reference cell coding, to standard 
#' reference cell coding.
#' 
#' @param est_xi Matrix of xi parameter estimates in mixture reference coding. KxQ
#' @return Returns vector `est_beta` of the probit regression coefficients
#' converted into reference cell coding with interactions. (K*Q)x1
#' 
#' @keywords internal
#' @export
convert_mix_to_ref <- function(est_xi) {
  est_beta <- est_xi
  for (i in 2:nrow(est_xi)) {
    est_beta[i, ] <- est_xi[i, ] - est_xi[1, ]
  }
  est_beta <- c(est_beta)
  return(est_beta)
  
  ## Alternative method
  # K <- nrow(est_xi)
  # Q <- ncol(est_xi)
  # est_beta <- numeric(K * Q)
  # # If only latent class as a covariate (no interactions)
  # if (Q == 1) {
  #   est_int_k1 <- NULL   # baseline class
  #   est_int_koth <- NULL # additional classes
  # # Position of interaction terms for additional covariates
  # } else if (Q > 1) {
  #   est_int_k1 <- K + 1:(Q-1)       # baseline class
  #   est_int_koth <- (K - 1) * (1:(Q-1)) # additional classes
  # }
  # 
  # # Get estimates for first class level, including all interactions
  # est_beta[c(1, est_int_k1)] <- est_xi[1, ]
  # # Get estimates for additional class levels, including all interactions
  # for (k in 2:K) {
  #   est_beta[c(k, (k + (Q-1)) + est_int_koth)] <- est_xi[k, ] - est_xi[1, ]
  # }
  # return(est_beta)
}

#' Convert from mixture reference coding to P(Y=1|-) conditional probabilities
#' 
#' @description
#' Convert regression estimates \eqn{\xi} from mixture reference coding to 
#' conditional probit regression probabilities, P(Y=1|-), for a given covariate.
#' 
#' @inheritParams swolca
#' @param est_xi Matrix of xi parameter estimates. KxQ
#' @param cov_name String vector of length 1 or 2 specifying the key covariate(s) 
#' for which to obtain the outcome probabilities. All covariate names must be 
#' included in `glm_form` and `V_data`. 
#' 
#' @return Returns dataframe `Phi_df` of the converted outcome probabilities 
#' (i.e., \eqn{\Phi} values) for the categories of the covariate(s) specified in 
#' `cov_name`. Number of rows is equal to the number of category combinations of 
#' the key covariate(s). The first K columns include the \eqn{\Phi} values 
#' for the K latent classes, evaluated at the covariate values listed in the 
#' remaining Q columns. The key covariate(s) are evaluated at all levels and 
#' the other covariates are evaluated at their reference levels. 
#' 
#' @importFrom stats terms as.formula model.matrix pnorm
#' @keywords internal
#' @export
#' 
convert_to_probs <- function(est_xi, glm_form, V_data, cov_name) {
  # check that cov_name is found in glm_form
  if (!all(sapply(cov_name, function(x) grepl(x, glm_form)))) {
    stop("all variables in cov_name must be specified in glm_form")
  } else if (!all(sapply(cov_name, function(x) x %in% colnames(V_data)))) {
    stop("all variables in cov_name must be found in V_data")
  } else if (grepl("c_all", glm_form)) {
    stop("glm_form must not contain the latent class variable c_all")
  }
  
  # Number of latent classes
  K <- nrow(est_xi)
  # Get all covariate names
  cov_names <- labels(stats::terms(stats::as.formula(glm_form)))
  # Get names of other covariates not highlighted in the plot
  oth_names <- cov_names[cov_names != cov_name]
  
  # Get levels for the key covariate(s)
  cov_levels <- lapply(cov_name, function(x) levels(V_data[[x]]))
  num_cov_levels <- nrow(expand.grid(cov_levels))
  # Get levels for the other covariates
  oth_levels <- lapply(oth_names, function(x) levels(V_data[[x]]))
  all_levels <- append(cov_levels, oth_levels)
  # All combinations of the levels of all covariates
  all_level_comb <- expand.grid(all_levels)
  colnames(all_level_comb) <- c(cov_name, oth_names)
  
  # Create design matrix from glm formula with new covariate level combinations
  all_model_mat <- stats::model.matrix(stats::as.formula(glm_form), all_level_comb)
  # Obtain Phi values corresponding to the new covariate level combinations
  all_Phi_df <- as.data.frame(sapply(1:K, function(k) 
    stats::pnorm(all_model_mat %*% est_xi[k, ])))
  colnames(all_Phi_df) <- paste0("Class", 1:K)
  # Rename the key covariate(s) for ease with plotting
  colnames(all_level_comb)[1:length(cov_name)] <- paste0("Cov", 1:length(cov_name))

  # Create new df of desired covariate combinations and corresponding Phi values
  Phi_df <- cbind(all_Phi_df, all_level_comb)
  # Restrict to levels of the key covariate(s) and key only the reference level 
  # values for the remaining covariates 
  Phi_df <- Phi_df[1:num_cov_levels, ]
  
  return(Phi_df)
  
}


#' Get confidence or credible interval
#' 
#' `get_ci` is a helper function that formats confidence or credible intervals
#' @param post_samp Numeric vector of posterior samples
#' @param lb Lower bound quantile of interval estimate. Default is 0.025 
#' corresponding to a 95% interval
#' @param ub Upper bound quantile of interval estimate. Default is 0.975 
#' corresponding to a 95% interval
#' @param digits Number of digits to round to. Default is 2.
#' @param string Boolean indicating if the interval should be returned as a 
#' string with parentheses. Default is `FALSE`, which returns the interval as a 
#' vector.
#' @importFrom stats quantile
#' @return Outputs string `ci` with formatted interval
#' @keywords internal
#' @export
get_ci <- function(post_samp, quant_lb = 0.025, quant_ub = 0.975, digits = 2, 
                   string = FALSE) {
  quantiles <- format(round(stats::quantile(post_samp, c(quant_lb, quant_ub)), 
                            digits), nsmall = digits)
  if (string) {
    ci <- paste0("(", quantiles[1], ", ", quantiles[2], ")")
  } else {
    ci <- as.numeric(quantiles)
  }
  return(ci)
}



#' Get posterior probability
#' 
#' `get_prob_pos` is a helper function that obtains the probability that the
#' posterior sample estimates are greater than some specified value
#' @param post_samp Numeric vector of posterior samples
#' @param cutoff Value specifying cutoff to compare posterior samples to. Default
#' is 0.
#' @param digits Number of digits to round to. Default is 2, where posterior 
#' probabilities smaller than 0.01 are denoted as "<0.01". 
#' @return Outputs `prob_pos` posterior probability
#' @keywords internal
#' @export
get_prob_pos <- function(post_samp, cutoff = 0, digits = 2) {
  prob_pos <- format(round(mean(post_samp > cutoff), digits), nsmall = digits)
  round_value <- 10^(-digits)
  if (prob_pos < round_value) {
    prob_pos <- paste0("<", round_value)
  }
  return(prob_pos)
}

