#===================================================
## Helper functions for WSOLCA, SOLCA, and WOLCA
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
#' (default), `solca`, or `wolca`
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
                         run_sampler = NULL, glm_form = NULL, K_max = NULL, 
                         adapt_seed = NULL, fixed_seed = NULL,
                         class_cutoff = NULL, alpha_adapt = NULL, 
                         eta_adapt = NULL, mu0_adapt = NULL, 
                         Sig0_adapt = NULL, K_fixed = NULL, alpha_fixed = NULL, 
                         eta_fixed = NULL, mu0_fixed = NULL, Sig0_fixed = NULL,
                         n_runs = NULL, burn = NULL, thin = NULL, 
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
    if (any(!is.null(c(alpha_adapt, eta_adapt)))) {
      if (length(alpha_adapt) != K_max) {
        stop("length of alpha_adapt must be the same as K_max")
      }
      if ((nrow(eta_adapt) != J) | (ncol(eta_adapt) != R)) {
        stop("eta_adapt must be a matrix with J rows and R columns, where J is
             the number of exposure items and R is the maximum number of 
             exposure categories")
      }
    }  
    if (any(!is.null(c(mu0_adapt, Sig0_adapt)))) {
      if (length(mu0_adapt) != K_max | !is.list(mu0_adapt) | 
          !(all(lapply(mu0_adapt, length) == q))) {
        stop("mu0_adapt must be a list of length K_max where each element is a 
           vector of length q (number of regression covariates excluding latent class)")
      }
      if (length(Sig0_adapt) != K_max | !is.list(Sig0_adapt) | 
          !(all(lapply(Sig0_adapt, nrow) == q)) | 
          !(all(lapply(Sig0_adapt, ncol) == q))) {
        stop("Sig0_adapt must be a list of length K_max where each element is a 
            qxq matrix, where q is the number of regression covariates excluding 
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
      if (any(!is.null(c(alpha_fixed, eta_fixed)))) {
        if (length(alpha_fixed) != K_fixed) {
          stop("length of alpha_fixed must be the same as K_fixed")
        }
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
      if (any(!is.null(c(mu0_fixed, Sig0_fixed)))) {
        if (length(mu0_fixed) != K_fixed | !is.list(mu0_fixed) | 
            !(all(lapply(mu0_fixed, length) == q))) {
          stop("mu0_fixed must be a list of length K_fixed where each element is a 
           vector of length q (number of regression covariates excluding latent class)")
        }
        if (length(Sig0_fixed) != K_fixed | !is.list(Sig0_fixed) | 
            !(all(lapply(Sig0_fixed, nrow) == q)) | 
            !(all(lapply(Sig0_fixed, ncol) == q))) {
          stop("Sig0_fixed must be a list of length K_fixed where each element is a 
            qxq matrix, where q is the number of regression covariates excluding 
           latent class)")
        }
      }
    }
    
    # Check MCMC parameters
    if (!all(is.null(c(n_runs, burn, thin)))) {
      if (!all(c(n_runs, burn, thin) %% 1 == 0) | !all(c(n_runs, burn, thin) >= 0)) {
        stop("n_runs, burn, and thin must be whole numbers")
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
#' @param q Number of regression covariates excluding class assignment
#' @param est_beta Vector of probit regression coefficients in reference cell 
#' coding. (K*q)x1. Order of betas must correspond to a formula with `c_all` as
#' the first covariate and all interaction terms involving `c_all` present.
#' @param ci_beta Matrix of interval estimates for probit regression 
#' coefficients in reference cell coding. (K*q)x2, where the first column is 
#' the lower bound and the second column is the upper bound. Set to `NULL` 
#' (default), if no interval estimate conversions are necessary.
#' 
#' @return Outputs list `xi_list` containing `est_xi`, a Kxq matrix of 
#' regression estimates in mixture reference coding. If `ci_beta` is not `NULL`,
#' `xi_list` also contains `est_xi_lb` and `est_xi_ub`, which are Kxq matrices
#' of the lower and upper bound interval estimates, respectively, in mixture 
#' reference coding.
#' 
#' @keywords internal
#' @seealso [fit_probit_wolca()]
#' @export
convert_ref_to_mix <- function(K, q, est_beta, ci_beta = NULL) {
  est_xi <- matrix(NA, nrow = K, ncol = q)
  if (!is.null(ci_beta)) {
    est_xi_lb <- est_xi_ub <- matrix(NA, nrow = K, ncol = q)
  }
  
  # If only latent class as a covariate (no interactions)
  if (q == 1) {
    est_int_k1 <- NULL   # baseline class
    est_int_koth <- NULL # additional classes
  # Position of interaction terms for additional covariates
  } else if (q > 1) {
    est_int_k1 <- K + 1:(q-1)       # baseline class
    est_int_koth <- (K - 1) * (1:(q-1)) # additional classes
  }
  
  # Get estimates for first class level, including all interactions
  est_xi[1, ] <- est_beta[c(1, est_int_k1)]
  if (!is.null(ci_beta)) {
    est_xi_lb[1, ] <- ci_beta[c(1, est_int_k1), 1]
    est_xi_ub[1, ] <- ci_beta[c(1, est_int_k1), 2]
  }
  # Get estimates for additional class levels, including all interactions
  for (k in 2:K) {
    est_xi[k, ] <- est_beta[c(k, (k + (q-1)) + est_int_koth)] + est_xi[1, ]
    if (!is.null(ci_beta)) {
      est_xi_lb[k, ] <- ci_beta[c(k, (k + (q-1)) + est_int_koth), 1] + est_xi_lb[1, ]
      est_xi_ub[k, ] <- ci_beta[c(k, (k + (q-1)) + est_int_koth), 2] + est_xi_ub[1, ]
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
#' @param est_xi Matrix of xi parameter estimates in mixture reference coding. Kxq
#' @return Returns vector `est_beta` of the probit regression coefficients
#' converted into reference cell coding with interactions. (K*q)x1
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
  # q <- ncol(est_xi)
  # est_beta <- numeric(K * q)
  # # If only latent class as a covariate (no interactions)
  # if (q == 1) {
  #   est_int_k1 <- NULL   # baseline class
  #   est_int_koth <- NULL # additional classes
  # # Position of interaction terms for additional covariates
  # } else if (q > 1) {
  #   est_int_k1 <- K + 1:(q-1)       # baseline class
  #   est_int_koth <- (K - 1) * (1:(q-1)) # additional classes
  # }
  # 
  # # Get estimates for first class level, including all interactions
  # est_beta[c(1, est_int_k1)] <- est_xi[1, ]
  # # Get estimates for additional class levels, including all interactions
  # for (k in 2:K) {
  #   est_beta[c(k, (k + (q-1)) + est_int_koth)] <- est_xi[k, ] - est_xi[1, ]
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
#' @param est_xi Matrix of xi parameter estimates. Kxq
#' @param cov_name String specifying name of covariate of interest. Must be 
#' included in `glm_form`.
#' @return Returns dataframe `probs` of the converted probabilities for the 
#' covariate specified in `cov_name`, with number of rows equal to K and number
#' of columns equal to the number of categories for the covariate (including the
#' baseline category) plus one. The first column specifies the latent class,
#' the second column corresponds to the baseline category intercept for all K 
#' latent classes, and the remaining columns correspond to the other categories 
#' for the covariate. 
#' 
#' @importFrom stats terms as.formula pnorm
#' @keywords internal
#' @export
#' 
convert_to_probs <- function(est_xi, glm_form, V_data, cov_name) {
  # check that cov_name is found in glm_form
  if (!grepl(cov_name, glm_form)) {
    stop("cov_name must be one of the variables specified in glm_form")
  } else if (grepl("c_all", glm_form)) {
    stop("glm_form must not contain the latent class variable c_all")
  }
  
  # Number of latent classes
  K <- nrow(est_xi)
  # Get covariate names
  cov_names <- labels(terms(as.formula(glm_form)))
  # Get index of covariate names corresponding to the covariate of interest
  select_cov <- which(cov_names == cov_name)
  
  # Get column indices for each variable in glm_form
  cov_col_inds <- attr(model.matrix(as.formula(glm_form), data = V_data), "assign")
  # Design matrix indices for covariate group, including intercept
  cols <- c(1, which(cov_col_inds == select_cov))
  
  # Get conversions for each category of the covariate group 
  probs <- as.data.frame(matrix(NA, nrow = K, ncol = (length(cols) + 1)))
  colnames(probs) <- c("Class", "Intercept", paste0(cov_name, 2:length(cols)))
  probs[, 1] <- 1:K
  for (categ in 1:length(cols)) {
    # Convert from mixture reference to probabilities
    probs[, categ + 1] <- stats::pnorm(est_xi[, 1] + 
                                         (categ > 1) * est_xi[, cols[categ]])
  }
  return(probs)
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
#' @importFrom stats quantile
#' @return Outputs string `ci` with formatted interval
#' @keywords internal
#' @export
get_ci <- function(post_samp, quant_lb = 0.025, quant_ub = 0.975, digits = 2) {
  quantiles <- format(round(stats::quantile(post_samp, c(quant_lb, quant_ub)), 
                            digits), nsmall = digits)
  ci <- paste0("(", quantiles[1], ", ", quantiles[2], ")")
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

