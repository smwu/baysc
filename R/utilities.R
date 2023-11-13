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


#' Convert from reference cell coding to factor reference coding 
#'
#' `convert_ref_to_comb` converts from reference cell coding to a combination of 
#' factor variable and reference cell coding, referred to as factor reference coding
#' 
#' @param beta_ref Matrix of probit coefficients in reference cell coding. (K*q)x1
#' @return Outputs Kxq matrix `beta_comb` of probit coefficients using factor 
#' variable coding
#' @keywords internal
#' @export
convert_ref_to_comb <- function(beta_ref) {
  beta_comb <- beta_ref
  for (i in 2:nrow(beta_ref)) {
    beta_comb[i, ] <- beta_ref[i, ] - beta_ref[1, ]
  }
  return(beta_comb)
}

#' Catch errors in input variables for external functions
#' 
#' @description
#' Catch input errors in variables necessary for package functions. All 
#' parameters are set to `NULL` by default so that error checks 
#' are only performed on relevant variables.
#' 
#' @inheritParams swolca
#' @return The function stops and an error message is displayed if the input 
#' variables are not acceptable
#' @details All parameters are set to `NULL` by default so that error checks 
#' are only performed on relevant variables.
#' 
#' @importFrom stats terms as.formula
#' @keywords internal
#' @export
catch_errors <- function(x_mat = NULL, y_all = NULL, sampling_wt = NULL, 
                         cluster_id = NULL, stratum_id = NULL, V = NULL,
                         run_sampler = NULL, glm_form = NULL, K_max = NULL, 
                         class_cutoff = NULL, alpha_adapt = NULL, 
                         eta_adapt = NULL, mu0_adapt = NULL, 
                         Sig0_adapt = NULL, K_fixed = NULL, alpha_fixed = NULL, 
                         eta_fixed = NULL, mu0_fixed = NULL, Sig0_fixed = NULL,
                         n_runs = NULL, burn = NULL, thin = NULL, 
                         save_res = NULL, save_path = NULL) {
  if (is.null(x_mat)) {
    stop("need to specify exposure matrix")
  } else {
    
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
      regr_vars <- labels(terms(as.formula(glm_form)))   
      if ("c_all" %in% regr_vars) {
        regr_vars <- regr_vars[regr_vars != "c_all"]
        if (!(all(regr_vars %in% colnames(V)))) {
          stop("all variables in glm_form except c_all must exist in V")
        }
      } else {
        if (!(all(regr_vars %in% colnames(V)))) {
          stop("all variables in glm_form must exist in V")
        }
      }
      if (substring(glm_form, 1, 1) != "~") {
        stop("glm_form must be a string starting with '~' that specifies a valid formula")
      }
    }
    
    # Check same number of individuals for x and y
    if (!is.null(y_all)) {
      if (n != length(y_all)) {
        stop("number of rows in x_mat must match length of y_all")
      }
    }
    
    # Check same number of individuals for x and sampling weights
    if (!is.null(sampling_wt)) {
      if (n != length(sampling_wt)) {
        stop("number of rows in x_mat must match length of sampling_wt")
      }
    }

    # If no clustering, assign each individual to their own cluster. Else, check
    # same number of individuals for x and clusters
    if (is.null(cluster_id)) {
      cluster_id <- 1:n
    } else if (n != length(cluster_id)) {
      stop("number of rows in x_mat must match length of cluster_id")
    }
    
    # Check same number of individuals for x and strata
    if (!is.null(stratum_id) & (n != length(stratum_id))) {
      stop("number of rows in x_mat must match length of stratum_id")
    }
    
    # Check same number of individuals for x and V
    if (!is.null(V)) {
      if (n != nrow(V)) {
        stop("number of rows in x_mat must match number of rows in V")
      }
      if (class(V)[1] != "data.frame") {
        stop("V must be a dataframe")
      }
    }
    
    # Check cutoff is between 0 and 1
    if (!is.null(class_cutoff)) {
      if (class_cutoff <= 0 | class_cutoff >= 1) {
        stop("class_cutoff must be a proportion in (0,1)")
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
    if (any(is.null(c(n_runs, burn, thin)))) {
      stop("n_runs, burn, and thin must be whole numbers")
    } else {
      if (!all(c(n_runs, burn, thin) %% 1 == 0) | !all(c(n_runs, burn, thin) >= 0)) {
        stop("n_runs, burn, and thin must be whole numbers")
      }
      if (burn > n_runs) {
        stop("n_runs must be larger than burn")
      }
    }

    # Check saving parameters
    if (!is.null(save_res)) {
      if (save_res) {
        if (is.null(save_path)) {
          stop("need to specify a path and file name in save_path, such as ~/Documents/run")
        }
      }
    }
  }
}


