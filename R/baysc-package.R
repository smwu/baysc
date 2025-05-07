#' @name baysc-package
#' @aliases baysc-package baysc
# #' @docType package
#' @title BAYesian Survey Clustering
#' 
#' @description
#' `baysc` is an R package for running Bayesian clustering methods on survey data.
#' A Bayesian latent class analysis (LCA) is available for eliciting underlying pattern 
#' profiles incorporating survey sampling weights and other survey design elements.
#'  
#' Options also exist for relating the pattern profiles to a binary outcome, either 
#' through a two-step approach after cluster is completed, or in a one-step approach 
#' where creation of the pattern profiles is directly informed by the outcome. 
#' Summary and plotting functions for visualizing output are available, as are 
#' diagnostic functions for examining convergence of the sampler. 
#' 
#' We provide an example dataset, `data_nhanes`, from the National Health and 
#' Nutrition Examination Survey (NHANES) that includes dietary intake and 
#' hypertension data for low-income women in the United States, as well as survey 
#' sampling weights and information on stratification and clustering in the 
#' survey design. A simulated dataset named `sim_data` is also provided to help 
#' with familiarization of the package.
#' 
#' @details
#' Use the [wolca()] function to run an unsupervised weighted overfitted LCA and 
#' obtain pattern profiles. [wolca_var_adjust()] provides a post-hoc variance 
#' adjustment that addresses variance underestimation. To examine the association 
#' of pattern profiles with a binary outcome through a two-step approach, run 
#' [wolca_svyglm()].
#' 
#' Use the [swolca()] function to run a supervised weighted overfitted LCA that 
#' allows information about the binary outcome to directly inform the creation 
#' of the pattern profiles. [swolca_var_adjust()] provides a post-hoc variance 
#' adjustment that addresses variance underestimation. 
#' 
#' More information about the methods can be found in the vignettes and in 
#' Wu et al. (2024).
#' 
#' @references 
#' Wu, S. M., Williams, M. R., Savitsky, T. D., & Stephenson, B. J. (2024). 
#' Derivation of outcome-dependent dietary patterns for low-income women 
#' obtained from survey data using a supervised weighted overfitted latent class 
#' analysis. Biometrics, 80(4), ujae122.
#' 
#' @keywords package
#' @importFrom RcppParallel RcppParallelLibs
#' 
#' @keywords internal
"_PACKAGE"