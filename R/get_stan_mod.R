#' Get SWOLCA_main stan model
#' 
#' Get stan model used in the `swolca()` function for the variance adjustment.
#' @param stan_name Name of stan model file without file extension
#' @return Returns `mod_stan` containing the specified stan model
#' @name rstantools_model_SWOLCA_main
#' @export
#' 
rstantools_model_SWOLCA_main <- function(stan_name) {
  mod_stan <- NULL
  if (stan_name == "SWOLCA_main") {
    mod_stan <- stanmodels$SWOLCA_main
  }
  return(mod_stan)
}

#' Get all stan models
#' 
#' Get all stan models available in the package.
#' @return Returns `stanmodels` object
#' @name stanmodels
#' @export
#' 
get_stanmodels <- function() {
  stanmodels
}

#' S4 class for SWOLCA_main stan model
#' 
#' S4 class for SWOLCA_main stan model inheriting `stanfit` class
#' @name Rcpp_rstantools_model_SWOLCA_main-class
#' @exportClass Rcpp_rstantools_model_SWOLCA_main
Rcpp_rstantools_model_SWOLCA_main <- setClass("Rcpp_rstantools_model_SWOLCA_main")

