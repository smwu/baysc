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

#' Obtains table of regression coefficients
#' 
#' @description
#' `get_regr_coefs` produces a summary table of the regression coefficients,
#' converted to standard reference cell coding. 
#' 
#' @inheritParams plot_theta_modes
#' @inheritParams plot_Phi_line
#' @param ci_level Numeric from 0 to 1 specifying the credible interval level. 
#' Default is 0.95, which gives a 95\% equal-tailed interval composed of the 
#' 2.5\% and 97.5\% quantiles. For `wolca()` results, this must match the 
#' `ci_level` parameter in the main function. 
#' @param digits Integer indicating the number of decimal places to be used. 
#' Default is 2, which rounds to the nearest hundredth. 
#' 
#' @return
#' Returns a character vector of the source code for regression coefficients 
#' table using the format specified in `format`.
#' 
#' @seealso [swolca()] [solca()] [wolca()] 
#' 
#' @importFrom dplyr mutate_if
#' @importFrom stats terms as.formula median quantile
#' @export
#'
#' #examples 
get_regr_coefs <- function(res, ci_level = 0.95, digits = 2) {
  
  if (!(ci_level > 0 & ci_level < 1)) {
    stop("ci_level must be between 0 and 1")
  }
  quant_lb <- (1 - ci_level) / 2
  quant_ub <- 1 - quant_lb
  
  # Estimates for `wolca()` can be obtained directly from the svyglm output
  if (is.null(res$estimates$xi_med)) {
    if (ci_level != res$data_vars$ci_level) {
      stop("ci_level must match the specified ci_level in the wolca() function")
    }
    beta <- as.data.frame(summary(res$estimates$fit)$coefficients)
    beta[, 4] <- ifelse(beta[, 4] < 10^(-digits), paste0("<", 10^(-digits)), 
                        format(round(beta[, 4], digits), digits))
    
  # Otherwise, for `swolca()` and `solca()`, create regression output after
  # converting from factor reference coding
  } else {
    # Obtain xi median and lower bound and upper bound estimates
    if (!is.null(res$estimates)) {
      # Estimates for `solca()` or adjusted estimates for `swolca()`
      est_xi <- res$estimates$xi_med
      est_lb <- apply(res$estimates$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_lb))
      est_ub <- apply(res$estimates$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_ub))
      est_red <- res$estimates$xi_red
    } else {
      # Unadjusted estimates for `swolca()`
      est_xi <- res$estimates_unadj$xi_med
      est_lb <- apply(res$estimates_unadj$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_lb))
      est_ub <- apply(res$estimates_unadj$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_ub))
      est_red <- res$estimates$xi_red
    }
    
    # Add outcome and latent class main and interaction terms to formula
    terms <- labels(stats::terms(stats::as.formula(res$data_vars$glm_form)))
    if (length(terms) > 0) {
      full_glm_form <- paste0("y_all ~ ", 
                              paste0("c_all * ", terms, collapse = " + ")) 
    } else {
      full_glm_form <- paste0("y_all ~ c_all") 
    }
    full_data <- data.frame(c_all = as.factor(res$estimates$c_all), 
                            res$data_vars$V_data)
    model_matrix <- model.matrix(as.formula(full_glm_form), data = full_data)
    
    K <- nrow(est_xi)
    q <- ncol(est_xi)
    beta <- as.data.frame(matrix(NA, nrow = ncol(model_matrix), ncol = 4))
    beta[, 1] <- colnames(model_matrix)
    colnames(beta) <- c("Covariate", "Estimate", "95% Cred Int",
                        "P(xi > 0)")
    
    # Intercept estimates
    beta[1, -1] <- c(est_xi[1, 1], get_ci(post_samp = est_red[, 1, 1]),
                     get_prob_pos(est_red[, 1, 1]))
    
    # Latent class main effect estimates
    for (i in 2:K) {
      beta[i, -1] <- c(stats::median(est_red[, i, 1] - est_red[, 1, 1]),
                       get_ci(est_red[, i, 1] - est_red[, 1, 1], digits = digits),
                       get_prob_pos(est_red[, i, 1] - est_red[, 1, 1], digits = digits))
    }
    
    # Additional covariates main effect estimates
    for (i in 2:q) {
      beta[K + (i-1), -1] <- c(est_xi[1, i], get_ci(est_red[, 1, i]),
                               get_prob_pos(est_red[, 1, i]))
    }
    
    # Additional covariates latent class interaction terms
    for (i in 2:q) {
      for (j in 2:K) {
        beta[q + (i-1)*(K-1) + (j-1), -1] <- 
          c(stats::median(est_red[, j, i] - est_red[, 1, i]),
            get_ci(est_red[, j, i] - est_red[, 1, i], digits = digits),
            get_prob_pos(est_red[, j, i] - est_red[, 1, i], digits = digits))
      }
    }
    beta$Estimate <- as.numeric(beta$Estimate)
  }
  
  # Print output
  beta <- dplyr::mutate_if(beta, is.numeric, round, digits = digits)
  beta
}