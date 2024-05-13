#' Reorder classes
#' 
#' @description
#' `reorder_classes` changes the order of the latent classes to help with 
#' plotting and to allow visualization of a different reference class for the 
#' regression coefficients.
#' 
#' @inheritParams plot_pattern_profiles
#' @param new_order Numeric vector specifying the new ordering of the classes. 
#' For example, if there are three classes and the new desired order is to have 
#' class 2 as the reference, followed by class 3 and class 1, then `new_order` 
#' would be `c(2, 3, 1)`.
#' 
#' @details
#' Latent class assignment variable `c_all` is also changed, so that if 
#' `new_order = c(2, 3, 1)`, then `c_all == 2` will become `c_all == 1`, 
#' `c_all = 3` will become `c_all = 2`, and `c_all == 1` will become `c_all = 3`.
#' 
#' @return
#' Returns object `res_new` that is identical to input `res` but has updated 
#' latent class ordering for `pi_red`, `theta_red`, `pi_med`, `theta_med`, and 
#' `c_all`. If `res` is a `"swolca"` object, `res_new` also includes updated 
#' latent class ordering for `xi_red` and `xi_med`. If `res` is a `"wolca"` 
#' object, reordering should be done prior to running [wolca_svyglm()] to obtain 
#' regression estimates with a different reference class. 
#' 
#' @seealso [plot_outcome_probs()] [plot_class_dist()] [summarize_res()]
#' 
#' @export
#'
#' @examples
#' # Load NHANES data
#' data(run_nhanes_swolca_results)
#' # Reorder latent classes
#' res_new <- reorder_classes(res = run_nhanes_swolca_results, 
#'                            new_order = c(3, 2, 5, 4, 1))
#' # Get posterior estimates for xi with class 3 as the reference
#' get_regr_coefs(res = res_new, ci_level = 0.95, digits = 2)                           
#' 
reorder_classes <- function(res, new_order) {
  # Check object class
  if (!inherits(res, c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  } else if ((inherits(res, "wolca")) & !is.null(res$estimates_svyglm)) {
    warning(paste0("For WOLCA, reordering of classes should be done before ",
                "calling wolca_svyglm(). res$estimates_svyglm should be NULL ",
                "prior to running this function."))
  }
  
  # Initialize res_new object
  res_new <- res
  if (!is.null(res$estimates_adjust)) {
    # Adjusted estimates 
    # Reorder classes for all pi and theta estimates
    res_new$estimates_adjust$pi_red <- res$estimates_adjust$pi_red[, new_order]
    res_new$estimates_adjust$theta_red <- res$estimates_adjust$theta_red[, , new_order, ]
    res_new$estimates_adjust$pi_med <- res$estimates_adjust$pi_med[new_order]
    res_new$estimates_adjust$theta_med <- res$estimates_adjust$theta_med[, new_order, ]
    # Reorder latent class assignment
    for (i in 1:5) {
      res_new$estimates_adjust$c_all[res$estimates_adjust$c_all == new_order[i]] <- i
    }
    if (is(res, "swolca")) {
      # If `swolca`, reorder classes for all xi estimates.
      res_new$estimates_adjust$xi_red <- res$estimates_adjust$xi_red[, new_order, ]
      res_new$estimates_adjust$xi_med <- res$estimates_adjust$xi_med[new_order, ]
    }
  } else {
    # Unadjusted estimates 
    # Reorder classes for all pi and theta estimates
    res_new$estimates$pi_red <- res$estimates$pi_red[, new_order]
    res_new$estimates$theta_red <- res$estimates$theta_red[, , new_order, ]
    res_new$estimates$pi_med <- res$estimates$pi_med[new_order]
    res_new$estimates$theta_med <- res$estimates$theta_med[, new_order, ]
    # Reorder latent class assignment
    for (i in 1:5) {
      res_new$estimates$c_all[res$estimates$c_all == new_order[i]] <- i
    }
    if (is(res, "swolca")) {
      # If `swolca`, reorder classes for all xi estimates.
      res_new$estimates$xi_red <- res$estimates$xi_red[, new_order, ]
      res_new$estimates$xi_med <- res$estimates$xi_med[new_order, ]
    }
  }
  
  return(res_new)
}


#' Obtains table of regression coefficients
#' 
#' @description
#' `get_regr_coefs` produces a summary table of the regression coefficients,
#' converted to standard reference cell coding. 
#' 
#' @inheritParams summarize_res
#' 
#' @details 
#' If `res` is a `swolca` object, any latent class can be chosen as the reference 
#' class level for which to display regression coefficients. Simply run 
#' [reorder_classes()] with the desired reference level as the first class in 
#' `new_order`, and then run [get_regr_coefs()] using the reordered `res` object.
#' 
#' If `res` is a `wolca` object, choosing a different reference level can only 
#' be done by re-running [wolca_svyglm()] with the new ordering. To do this, run 
#' [reorder_classes()], then run [wolca_svyglm()] to obtain regression 
#' estimates with a different reference class, and then finally use the 
#' [get_regr_coefs()] function. 
#' 
#' @return
#' Returns dataframe `beta` containing the following columns:
#' \describe{
#'   \item{\code{Covariate}}{Names of covariate terms.}
#'   \item{\code{Estimate}}{Regression estimates in reference cell coding format.}
#'   \item{\code{LB}}{Regression estimate lower bound corresponding to the level 
#'   specified in `ci_level`. Calculated using posterior samples for `swolca` 
#'   objects and confidence intervals for `wolca` objects.}
#'   \item{\code{UB}}{Regression estimate upper bound corresponding to the level 
#'   specified in `ci_level`. Calculated using posterior samples for `swolca` 
#'   objects and confidence intervals for `wolca` objects.}
#'   \item{\code{P(xi>0)} or \code{p-value}}{For `swolca` objects, the proportion 
#'   of posterior samples greater than 0. For `wolca` objects, the p-value 
#'   obtained from the `wolca_svyglm()` function}
#' }
#' 
#' @seealso [reorder_classes()] [plot_outcome_probs()] [summarize_res()]
#' 
#' @importFrom dplyr mutate_if
#' @importFrom stats terms as.formula median quantile
#' @export
#'
#' @examples
#' data(run_nhanes_swolca_results)
#' get_regr_coefs(res = run_nhanes_swolca_results, ci_level = 0.95, digits = 2)
#' 
get_regr_coefs <- function(res, ci_level = 0.95, digits = 2) {
  
  # Check object class
  if (!(inherits(res, c("swolca", "wolca")))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  if (!(ci_level > 0 & ci_level < 1)) {
    stop("ci_level must be between 0 and 1")
  }
  quant_lb <- (1 - ci_level) / 2
  quant_ub <- 1 - quant_lb
  
  # Set pointer to adjusted or unadjusted estimates
  if (!is.null(res$estimates_adjust)) {
    # Adjusted estimates
    estimates <- res$estimates_adjust
  } else {
    # Unadjusted estimates
    estimates <- res$estimates
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
                          y_all = res$data_vars$y_all,
                          res$data_vars$V_data)
  model_matrix <- model.matrix(as.formula(full_glm_form), data = full_data)
  beta <- as.data.frame(matrix(NA, nrow = ncol(model_matrix), ncol = 5))
  beta[, 1] <- colnames(model_matrix)
  
  # Estimates for `wolca()` can be obtained directly from the svyglm output
  if (is(res, "wolca")) {
    if (ci_level != res$data_vars$ci_level) {
      stop("ci_level must match the specified ci_level in the wolca() function")
    }
    colnames(beta) <- c("Covariate", "Estimate", "LB", "UB", "p-value")
    beta[, c(2, 5)] <- res$estimates_svyglm$fit_summary$coefficients[, c(1, 4)]
    beta[, 2] <- format(round(beta[, 2], digits), digits)  
    beta[, 3] <- format(round(convert_mix_to_ref(res$estimates_svyglm$xi_est_lb), 
                              digits), digits)
    beta[, 4] <- format(round(convert_mix_to_ref(res$estimates_svyglm$xi_est_ub),
                              digits), digits)
    beta[, 5] <- ifelse(beta[, 5] < 10^(-digits), paste0("<", 10^(-digits)), 
                        format(round(beta[, 5], digits), digits))
    
  # Otherwise, for `swolca()`, create regression output after converting from 
  # mixture reference coding
  } else {
    # Obtain xi median and lower bound and upper bound estimates
    est_xi <- estimates$xi_med
    est_lb <- apply(estimates$xi_red, c(2, 3),
                    function(x) stats::quantile(x, quant_lb))
    est_ub <- apply(estimates$xi_red, c(2, 3),
                    function(x) stats::quantile(x, quant_ub))
    est_red <- estimates$xi_red

    K <- nrow(est_xi)
    Q <- ncol(est_xi)
    colnames(beta) <- c("Covariate", "Estimate", "LB", "UB", "P(xi > 0)")
    
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
    for (i in 2:Q) {
      beta[K + (i-1), -1] <- c(est_xi[1, i], get_ci(est_red[, 1, i]),
                               get_prob_pos(est_red[, 1, i]))
    }
    
    # Additional covariates latent class interaction terms
    for (i in 2:Q) {
      for (j in 2:K) {
        beta[Q + (i-1)*(K-1) + (j-1), -1] <- 
          c(stats::median(est_red[, j, i] - est_red[, 1, i]),
            get_ci(est_red[, j, i] - est_red[, 1, i], digits = digits),
            get_prob_pos(est_red[, j, i] - est_red[, 1, i], digits = digits))
      }
    }
    beta$Estimate <- as.numeric(beta$Estimate)
    beta$LB <- as.numeric(beta$LB)
    beta$UB <- as.numeric(beta$UB)
  }
  
  # Print output
  beta <- dplyr::mutate_if(beta, is.numeric, round, digits = digits)
  beta
}





#' Get levels of a variable 
#'
#' `get_levels` obtains the levels of a specified variable in a dataframe. For 
#' factor varialbes, the levels are directly returned. For continuous variables,
#' the levels are returned as an empty string. Internal function used in 
#' [vars_across_class()].
#' 
#' @param df Dataframe including variable of interest
#' @param var String specifying variable of interest
#' 
#' @return Outputs string vector `levels` containing the levels of the variable 
#' of interest. If the variable is not a factor, outputs the empty string `""`.
#' 
#' @keywords internal
#' @seealso [vars_across_class()]
#' @export
get_levels <- function(df, var) {
  # For factor variables, get levels directly
  if (is.factor(df[[var]])) {
    levels <- levels(df[[var]])
    # For continuous variables, set levels to empty string
  } else {
    levels <- ""
  }
  return(levels)
}

#' Get the covariate distribution across latent classes
#' 
#' @description 
#' `get_cov_props` gives the column means or percentages across latent classes 
#' for a given demographic variable, calculated for the population using survey 
#' weights and methods from the `survey` R package. Internal function used in 
#' [vars_across_class()].
#' 
#' @inheritParams vars_across_class
#' @param svy_design Survey design defined using the `svydesign()` function in 
#' the `survey` package
#' @param cov String specifying covariate of interest. Use `"population"` or 
#' `"sample"` to get the class-specific population or sample sizes, respectively.
#' @param var_levels String vector specifying levels of the covariate of interest. 
#' If the covariate is continuous, set this to the empty string `""`.
#' @param res Results from running [swolca()] or [wolca()], with or without the 
#' variance adjustment. This is necessary when `cov == "population"` but can 
#' otherwise be set to the default value of `NULL`.
#' 
#' @details
#' If the covariate is a factor, the function outputs percentages for each level. 
#' If the covariate is continuous, the function outputs the mean for each level. 
#' If the covariate is the population class membership proportions (i.e., 
#' `cov == "population"`), the posterior standard deviation is also provided 
#' using the results from running [swolca()] or [wolca()]. The standard errors 
#' for the other variables are not provided because they would be underestimates 
#' due to not accounting for variability in the latent class assignments. If the 
#' covariate is the sample class membership proportions (i.e., `cov == "sample"`), 
#' then the percentages provided are exact.
#' 
#' @return
#' Outputs `output` dataframe with number of rows equal to the number of levels 
#' in the covariate of interest (one row if the covariate is continuous), and 
#' number of columns equal to the number of latent classes, plus an additional 
#' column for the overall value over all classes. 
#' 
#' @importFrom survey svydesign svymean svyby
#' @importFrom stats as.formula sd
#' @keywords internal
#' @export
get_cov_props <- function(svy_design, cov, var_levels, col_props = TRUE, 
                          res = NULL, digits = 1) {
  # Check object class
  if (!is.null(res)) {
    if (!(inherits(res, c("swolca", "wolca")))) {
      stop("res must be NULL or an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
    }
  }
  # Number of covariate categories
  n_levels <- length(var_levels)
  # Number of columns is number of latent classes plus one for "Overall"
  n_cols <- length(table(svy_design$variables$Class)) + 1
  # Initialize output table
  output <- as.data.frame(matrix(NA, nrow = n_levels, ncol = n_cols))
  colnames(output) <- c(paste0("Pattern", 1:(n_cols-1)), "Total")
  
  # If sample n, just tabulate number of individuals in each latent class. No SE
  if (cov == "sample") {
    temp_tot <- length(svy_design$variables$Class)
    temp_class <- table(svy_design$variables$Class)
    output[, n_cols] <- ""
    output[, -n_cols] <- format(round(temp_class / temp_tot * 100, digits), 
                                nsmall = digits)
  # If population N, estimate pop mean and SD for each class using posterior pi
  } else if (cov == "population") {
    output[, n_cols] <- ""
    # Get posterior SD of the class proportions, using adjusted or unadjusted
    if (!is.null(res$estimates_adjust)) {
      output[, -n_cols] <- unlist(apply(res$estimates_adjust$pi_red, 2, function(x) 
        paste0(format(round(mean(x)*100, digits), nsmall = digits), " (", 
               format(round(stats::sd(x)*100, digits), nsmall = digits), ")")))
    } else {
      output[, -n_cols] <- unlist(apply(res$estimates$pi_red, 2, function(x) 
        paste0(format(round(mean(x)*100, digits), nsmall = digits), " (", 
               format(round(stats::sd(x)*100, digits), nsmall = digits), ")")))
    }
  # If continuous covariate, estimate population mean for each class
  } else if (!is.factor(svy_design$variables[[cov]])) {
    temp_tot <- as.data.frame(survey::svymean(stats::as.formula(paste0("~", cov)), 
                                              svy_design))
    temp_class <- as.data.frame(survey::svyby(stats::as.formula(paste0("~", cov)), 
                                              ~Class, svy_design, survey::svymean))
    output[, n_cols] <- format(round(temp_tot$mean, digits), nsmall = digits)
    output[, -n_cols] <- sapply(1:(n_cols-1), function(x)
      format(round(temp_class[x, 1 + 1:n_levels], digits), nsmall = digits))
  # If factor covariate, estimate population column or row proportions
  } else {
    if (col_props) {
      # column proportions (category proportions for a given class)
      temp_tot <- as.data.frame(survey::svymean(stats::as.formula(paste0("~", cov)), 
                                                svy_design))
      temp_class <- as.data.frame(survey::svyby(stats::as.formula(paste0("~", cov)), 
                                                ~Class, svy_design, survey::svymean))
      output[, n_cols] <- format(round(temp_tot$mean*100, digits), nsmall = digits)
      output[, -n_cols] <- as.data.frame(sapply(1:(n_cols-1), function(x)
        format(round(temp_class[x, 1 + 1:n_levels]*100, digits), nsmall = digits)))
    } else {
      # row proportions (class proportions for a given category)
      temp_class <- as.data.frame(survey::svyby(~Class,  stats::as.formula(paste0("~", cov)), 
                                                svy_design, survey::svymean))
      output[, n_cols] <- 100
      output[, -n_cols] <- format(round(temp_class[1:n_levels, 2:n_cols]*100, 
                                       digits), nsmall = digits)
    }
  }
  return(output)
}


#' Create table of class distributions across variables
#' 
#' @description
#' `vars_across_class` displays the distribution of variables across latent 
#' classes in the population. 
#' 
#' @inheritParams swolca
#' @param c_all nx1 factor vector of latent class assignments, typically obtained 
#' from the [swolca()], [wolca()], [swolca_var_adjust()] or [wolca_var_adjust()]
#' functions. Factor levels should be labeled with the names that are to appear 
#' in the output table.
#' @param cov_df Dataframe with n rows, with columns consisting of the variables 
#' that are to be shown with their distributions across latent classes. Factor 
#' variables should be labeled with the names and levels that are to appear in the 
#' output table. 
#' @param col_props Boolean indicating if factor variables should have percentages 
#' reported as column totals (default `TRUE`) that provide the percentage of the 
#' population in each category for a given latent class, or reported as row 
#' totals (`FALSE`) that provide the percentage of the population in each latent 
#' class for a given category. 
#' @param digits Integer indicating the number of decimal places to be used. 
#' Default is 1, which rounds to the nearest tenth. 
#' @param res Results from [swolca()], [wolca()], [swolca_var_adjust()] or 
#' [wolca_var_adjust()].
#' 
#' @details
#' The dataframe output includes a first row that presents the population 
#' percentage of individuals falling into each latent class, as well as the 
#' corresponding posterior standard deviation, obtained from the \eqn{\pi} 
#' parameters from a `"swolca"` or `"wolca"` object. The second row is the sample 
#' percentage falling into each latent class. 
#' 
#' Subsequent rows for continuous variables display the mean values for each 
#' latent class in the population, adjusting for survey design. Subsequent rows 
#' for factor variables display either column totals (percentage of the 
#' population in each category for a given latent class) or row totals 
#' (percentage of the population in each latent class for a given category), 
#' depending on the `col_props` input variable, adjusting for survey design. 
#' Survey design adjustments are carried out post-hoc using methods provided in 
#' the `survey` R package (Lumley, 2004). Since these methods do not account for 
#' variability in the latent classes, their standard error estimates will be 
#' underestimates and are not reported. 
#' 
#' @references
#' Lumley T (2004). “Analysis of Complex Survey Samples.” Journal of Statistical 
#' Software, 9(1), 1-19. R package verson 2.2.
#' 
#' @return 
#' Outputs dataframe `output_df` displaying the distribution of variables across 
#' the latent classes in the population. `output_df` has columns for the 
#' variable names, the category names for factor variables, the distribution of 
#' the variables across the latent classes, and the distribution of the 
#' variables overall in the population.
#' 
#' @export
#' @examples
#' data("data_nhanes")
#' data("run_nhanes_swolca_results")
#' res <- run_nhanes_swolca_results
#' c_all <- factor(res$estimates_adjust$c_all, levels = 1:5, 
#'                 labels = paste0("C", 1:5))
#' cov_df <- dplyr::select(data_nhanes, RIDAGEYR, age_cat, racethnic, smoker, 
#'                         physactive)
#' cov_df <- dplyr::rename(cov_df, Age = RIDAGEYR)
#' cov_df <- dplyr::mutate(cov_df, 
#'    Age_Group = factor(age_cat, levels = c(1, 2, 3), 
#'                       labels = c("[20, 40)", "[40, 60)", ">=60")),
#'    Race_Ethnicity = factor(racethnic, c(1, 2, 3, 4, 5),
#'                       labels = c("NH White", "NH Black", "NH Asian", 
#'                                  "Hispanic/Latino", "Other/Mixed")),
#'    Current_Smoker = factor(smoker, levels = c(0, 1), labels = c("No", "Yes")),
#'    Physical_Activity = factor(physactive, levels = c("Inactive", "Active")),
#'    .keep = "unused")
#' output_df <- vars_across_class(c_all = c_all, cov_df = cov_df, 
#'                               sampling_wt = data_nhanes$sample_wt, 
#'                               stratum_id = data_nhanes$stratum_id,
#'                               cluster_id = data_nhanes$cluster_id,
#'                               digits = 1, col_props = TRUE, res = res)
#'                  
vars_across_class <- function(c_all, cov_df, sampling_wt, stratum_id, cluster_id, 
                             digits = 1, col_props = TRUE, res) {
  if (!is.factor(c_all)) {
    stop("c_all must be a factor")
  }
  # Check object class
  if (!(inherits(res, c("swolca", "wolca")))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  # Covariate variable names
  var_names <- colnames(cov_df)
  # Covariate variable levels
  var_levels_list <- sapply(var_names, function(x) get_levels(df = cov_df, var = x))
  # Continuous variables
  cont_vars <- which(sapply(cov_df, function(x) !is.factor(x)))
  # Get number of rows and columns for the dataframe
  num_rows <- length(unlist(sapply(cov_df, levels))) + length(cont_vars) + 2
  num_cols <- length(levels(c_all)) + 3
  
  # Set survey design
  svy_design <- survey::svydesign(id = ~cluster_id,
                                  weights = ~sampling_wt,
                                  strata=~stratum_id,
                                  data = data.frame(cov_df, Class = c_all))
  
  # Initialize dataframe
  output_df <- as.data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  # Column names, including latent class names
  colnames(output_df) <- c("Variable", "Level", levels(c_all), "Overall")
  # Variable names
  output_df[, 1] <- c("N: % (posterior SD %)", "n: %", 
     unlist(sapply(1:length(var_names), function(x) 
        if(x %in% cont_vars){  # continuous variables
          c(paste0(var_names[x], ": Mean"))
        } else {  # factor variables
          c(paste0(var_names[x], ": %"), rep("", length(var_levels_list[[x]]) - 1))
        })))
  # Variable levels
  output_df[, 2] <- c("", "", unlist(var_levels_list))
  # Get estimates for population and sample sizes
  output_df[1, -c(1, 2)] <- unlist(get_cov_props(svy_design = svy_design, 
                                         cov = "population", var_levels = "", 
                                         res = res, digits = digits))
  output_df[2, -c(1, 2)] <- unlist(get_cov_props(svy_design = svy_design, 
                                                 cov = "sample", var_levels = "", 
                                                 digits = digits))
  # Get estimates for covariates
  row <- 2
  for (i in 1:length(var_names)) {
    var_levels <- var_levels_list[[i]]
    output_df[row + 1:length(var_levels), -c(1, 2)] <- 
      unlist(get_cov_props(svy_design = svy_design, cov = var_names[i], 
                            var_levels = var_levels, digits = digits, 
                            col_props = col_props))
    row <- row + length(var_levels)
  }
  # Output table
  return(output_df)
}


#' Obtains table of regression coefficients
#' 
#' @description
#' `summarize_res` produces a summary table of the regression coefficients,
#' converted to standard reference cell coding. 
#' 
#' @inheritParams plot_pattern_profiles
#' @param ci_level Numeric from 0 to 1 specifying the credible interval level. 
#' Default is 0.95, which gives a 95\% equal-tailed interval composed of the 
#' 2.5\% and 97.5\% quantiles. For `wolca` objects, this must match the 
#' `ci_level` parameter in the main function. 
#' @param digits Integer indicating the number of decimal places to be used. 
#' Default is 2, which rounds to the nearest hundredth. 
#' 
#' @return
#' Returns dataframe `estimates_df` of the coefficients estimates and their 
#' lower and upper bounds using the `ci_level` specified, rounded to the number 
#' of digits specified in `digits`. \eqn{\pi} (class membership probability) 
#' parameters are listed first, followed by \eqn{\theta} (item level probability) 
#' parameters and \eqn{\xi} (regression coefficient) parameters. 
#' 
#' @seealso [get_regr_coefs()] 
#' 
#' @importFrom dplyr mutate_if
#' @importFrom stats terms as.formula median quantile
#' @export
#'
#' @examples
#' data(run_nhanes_swolca_results)
#' estimates_df <- summarize_res(res = run_nhanes_swolca_results, ci_level = 0.95, 
#'                               digits = 2)
#' 
summarize_res <- function(res, ci_level = 0.95, digits = 2) {
  
  # Check object class
  if (!(inherits(res, c("swolca", "wolca")))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  if (!(ci_level > 0 & ci_level < 1)) {
    stop("ci_level must be between 0 and 1")
  }
  quant_lb <- (1 - ci_level) / 2
  quant_ub <- 1 - quant_lb
  
  # Set pointer to adjusted or unadjusted estimates
  if (!is.null(res$estimates_adjust)) {
    # Adjusted estimates
    estimates <- res$estimates_adjust
  } else {
    # Unadjusted estimates
    estimates <- res$estimates
  }
  
  # Check WOLCA for ci_level error and obtain Q
  if (is(res, "wolca")) {  # `wolca()`
    if (ci_level != res$data_vars$ci_level) {
      stop("ci_level must match the specified ci_level in the wolca() function")
    }
    # Obtain number of regression parameters excluding latent class
    Q <- dim(estimates$xi_est)[2]
  } else {  # `swolca()
    # Obtain number of regression parameters excluding latent class
    Q <- dim(estimates$xi_med)[2]
  }
  
  # Obtain dimensions
  K <- length(estimates$pi_med)
  J <- dim(estimates$theta_med)[1]
  R <- dim(estimates$theta_med)[3]
  
  # Obtain median and lower bound and upper bound estimates
  # pi
  pi_med <- estimates$pi_med
  pi_lb <- apply(estimates$pi_red, 2,
                 function(x) stats::quantile(x, quant_lb))
  pi_ub <- apply(estimates$pi_red, 2,
                 function(x) stats::quantile(x, quant_ub))
  pi_summary <- cbind(pi_med, pi_lb, pi_ub)
  # theta
  theta_med <- c(estimates$theta_med) # iterates over j, then k, then r
  theta_lb <- c(apply(estimates$theta_red, c(2, 3, 4), 
                      function(x) stats::quantile(x, quant_lb)))
  theta_ub <- c(apply(estimates$theta_red, c(2, 3, 4), 
                      function(x) stats::quantile(x, quant_ub)))
  theta_summary <- cbind(theta_med, theta_lb, theta_ub)
  # xi
  # Estimates for `wolca()` can be obtained directly from output
  if (is(res, "wolca")) {
    xi_med <- c(res$estimates$xi_est)
    xi_lb <- c(res$estimates$xi_est_lb)
    xi_ub <- c(res$estimates$xi_est_ub)
    xi_summary <- cbind(xi_med, xi_lb, xi_ub)
  } else {  # `swolca()` 
    xi_med <- c(estimates$xi_med)  # iterates over k then q
    xi_lb <- c(apply(estimates$xi_red, c(2, 3),
                     function(x) stats::quantile(x, quant_lb)))
    xi_ub <- c(apply(estimates$xi_red, c(2, 3),
                     function(x) stats::quantile(x, quant_ub)))
    xi_summary <- cbind(xi_med, xi_lb, xi_ub)
  }
  
  num_est <- sum(c(length(pi_med), length(theta_med), length(xi_med)))
  estimates_df <- as.data.frame(matrix(NA, nrow = num_est, ncol = 4))
  colnames(estimates_df) <- c("Parameter", "Estimate", "Lower Bound", "Upper Bound")
  # Get vector list of names for each of the parameters
  pi_names <- paste0("pi_", 1:K)
  # For theta, iterates over j, then k, then r
  theta_names <- paste0("theta_", 
                        paste0(paste0(rep(1:J, times = K*R), "_", 
                                      rep(1:K, each = J, times = R)), "_", 
                               rep(1:R, each = J*K))) 
  # For xi, iterates over k, then r
  xi_names <- paste0("xi_", paste0(rep(1:K, times = Q), "_",
                                   rep(1:Q, each = K)))
  estimates_df[, 1] <- c(pi_names, theta_names, xi_names)
  estimates_df[, -1] <- round(rbind(pi_summary, theta_summary, xi_summary), 
                              digits = digits)
  
  # Return dataframe of estimates
  return(estimates_df)
}


#' Obtain parameter MCMC iterations
#' 
#' @description
#' `get_param_mcmc` extracts the MCMC iteration values for all parameters 
#' obtained through the MCMC sampler.
#' 
#' @inheritParams plot_pattern_profiles
#' 
#' @return
#' Returns list `param_mcmc` containing the following:
#' \describe{
#'   \item{\code{pi_mcmc}}{MxK dataframe of the class membership probability 
#'   parameters, \eqn{\pi}, where M is the number of iterations and K is the 
#'   number of latent classes.}
#'   \item{\code{theta_mcmc}}{Mx(JxKxR) dataframe of the item level probability 
#'   parameters, \eqn{\theta}, where J is the number of items and R is the 
#'   maximum number of item levels.}
#' }
#' If `res` is a `"swolca"` object, `param_mcmc` also contains `xi_mcmc`, a 
#' Mx(KxQ) dataframe of the regression parameters, \eqn{\xi}, where Q is the 
#' number of covariates, excluding latent class indicators, in the regression. 
#' 
#' @seealso [summarize_res()] 
#' @export
#'
#' @examples
#' data(run_nhanes_swolca_results)
#' param_mcmc <- get_param_mcmc(res = run_nhanes_swolca_results)
#' 
get_param_mcmc <- function(res) {
  # Check object class
  if (!(inherits(res, c("swolca", "wolca")))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  # Set pointer to adjusted or unadjusted estimates
  if (!is.null(res$estimates_adjust)) {
    # Adjusted estimates
    estimates <- res$estimates_adjust
  } else {
    # Unadjusted estimates
    estimates <- res$estimates
  }
  
  # Obtain dimensions
  K <- length(estimates$pi_med)
  J <- dim(estimates$theta_med)[1]
  R <- dim(estimates$theta_med)[3]
  
  # Get vector list of names for each of the parameters
  pi_names <- paste0("pi_", 1:K)
  # For theta, iterates over j, then k, then r
  theta_names <- paste0("theta_", 
                        paste0(paste0(rep(1:J, times = K*R), "_", 
                                      rep(1:K, each = J, times = R)), "_", 
                               rep(1:R, each = J*K))) 
  if (is(res, "swolca")) {
    Q <- dim(estimates$xi_med)[2]
    # For xi, iterates over k, then r
    xi_names <- paste0("xi_", paste0(rep(1:K, times = Q), "_",
                                     rep(1:Q, each = K)))
  }
  
  # Get MCMC iterations for all estimates
  pi_mcmc <- as.data.frame(estimates$pi_red)
  colnames(pi_mcmc) <- pi_names
  theta_mcmc <- as.data.frame(estimates$theta_red)
  colnames(theta_mcmc) <- theta_names
  param_mcmc <- list(pi_mcmc = pi_mcmc, theta_mcmc = theta_mcmc)
  
  # Add xi if swolca
  if (is(res, "swolca")) {
    xi_mcmc <- as.data.frame(estimates$xi_red)
    colnames(xi_mcmc) <- xi_names
    param_mcmc$xi_mcmc = xi_mcmc
  }
  
  return(param_mcmc)
}


#' Obtain DIC-6
#' 
#' @description
#' `get_dic` returns the DIC-6 for the model to be used as a measure of model fit.
#' 
#' @inheritParams plot_pattern_profiles
#' 
#' @return
#' Returns the DIC-6.
#' 
#' @details
#' The deviance information criterion (DIC) can be used as a metric for model 
#' goodness of fit for mixture models (Spiegelhater et al., 2002; Celeux et al., 2006). 
#' It compares two values:
#' 1) the median of the likelihood over all iterations
#' 2) the likelihood calculated with the posterior median estimates of all parameters.
#' We use the DIC-6, which increases the penalty of complexity to reduce 
#' overfitting (Stephenson et al., 2022). The DIC-6 is calculated as 
#' \eqn{3\bar{D}(\theta)-2D(\bar{\theta})) 
#' = -6E[\log L(D|\theta)|D] + 4\log L(y|\bar{\theta})}.
#' where \eqn{E[\cdot]} denotes the expectation, \eqn{\bar{D}(\theta)} is the 
#' posterior median observed deviance and \eqn{D(\bar{\theta})} is the deviance 
#' of the posterior median.
#'  
#' @seealso [summarize_res()] 
#' @export
#' @references 
#' Celeux, G., Forbes, F., Robert, C. P., Titterington, D. M. et al. (2006) 
#' Deviance information criteria for missing data models. Bayesian Analysis, 1, 651-673.
#' 
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and Van Der Linde, A. (2002) 
#' Bayesian measures of model complexity and fit. Journal of the Royal 
#' Statistical Society: Series B (Statistical Methodology), 64, 583-639.
#' 
#' Stephenson, B. J., Herring, A. H., and Olshan, A. F. (2022). Derivation of 
#' maternal dietary patterns accounting for regional heterogeneity. Journal of 
#' the Royal Statistical Society Series C: Applied Statistics 71, 1957–1977.
#'
#' @examples
#' data(run_nhanes_swolca_results)
#' dic <- get_dic(res = run_nhanes_swolca_results)
#' 
get_dic <- function(res) {
  # Check object class
  if (!(inherits(res, c("swolca", "wolca")))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  # Set pointer to adjusted or unadjusted estimates
  if (!is.null(res$estimates_adjust)) {
    # Adjusted estimates
    estimates <- res$estimates_adjust
  } else {
    # Unadjusted estimates
    estimates <- res$estimates
  }
  
  dic <- -6 * median(res$MCMC_out$loglik_MCMC) + 4 * sum(estimates$log_lik_med)
  
  return(dic)
}






