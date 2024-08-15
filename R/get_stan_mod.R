#' Get SWOLCA_main stan model
#' 
#' Get stan model used in the `swolca()` function for the variance adjustment.
#' @return Returns the specified stan model
#' @name rstantools_model_SWOLCA_main
#' @keywords internal
#' @export
#' 
rstantools_model_SWOLCA_main <- function() {
  stanmodels$SWOLCA_main
}

#' Get WOLCA_main stan model
#' 
#' Get stan model used in the `wolca()` function for the variance adjustment.
#' @return Returns the specified stan model
#' @name rstantools_model_WOLCA_main
#' @keywords internal
#' @export
#' 
rstantools_model_WOLCA_main <- function() {
  stanmodels$WOLCA_main
}

#' Get all stan models
#' 
#' Get all stan models available in the package.
#' @return Returns `stanmodels` object
#' @name stanmodels
#' @importFrom rstantools rstan_config
#' @keywords internal
#' @export
#' 
get_stanmodels <- function() {
  stanmodels
}

#' S4 class for SWOLCA_main stan model
#' 
#' S4 class for SWOLCA_main stan model inheriting `stanfit` class
#' @name Rcpp_rstantools_model_SWOLCA_main-class
#' @keywords internal
#' @exportClass Rcpp_rstantools_model_SWOLCA_main
Rcpp_rstantools_model_SWOLCA_main <- setClass("Rcpp_rstantools_model_SWOLCA_main")

#' S4 class for WOLCA_main stan model
#' 
#' S4 class for WOLCA_main stan model inheriting `stanfit` class
#' @name Rcpp_rstantools_model_WOLCA_main-class
#' @keywords internal
#' @exportClass Rcpp_rstantools_model_WOLCA_main
Rcpp_rstantools_model_WOLCA_main <- setClass("Rcpp_rstantools_model_WOLCA_main")

