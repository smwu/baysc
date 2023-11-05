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

