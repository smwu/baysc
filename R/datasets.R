#' Simulated data
#' 
#' A dataset simulated using `simulate_data.R`
#' @docType data
#' @usage data(sim_data)
#' @format A list with 20 elements:
#' \describe{
#'   \item{\code{samp_ind}}{Vector of population indices for sampled individuals. nx1}
#'   \item{\code{sample_wt}}{Vector of individual survey sampling weights. nx1}
#'   \item{\code{true_Si}}{Vector of stratum indicators. nx1}
#'   \item{\code{cluster_id}}{Vector of cluster indicators. nx1}
#'   \item{\code{X_data}}{Matrix of multivariate categorical exposure. nxp}
#'   \item{\code{Y_data}}{Vector of outcomes. nx1}
#'   \item{\code{N}}{Population size}
#'   \item{\code{N_s}}{Vector of stratum-specific population sizes. Sx1}
#'   \item{\code{p}}{Number of exposure items}
#'   \item{\code{d}}{Number of exposure categories}
#'   \item{\code{S}}{Number of strata}
#'   \item{\code{true_K}}{True number of latent classes}
#'   \item{\code{true_pi_s}}{Matrix of true class membership probabilities by stratum. SxK}
#'   \item{\code{true_pi}}{Vector of true class membership probabilities across strata. Kx1}
#'   \item{\code{true_global_patterns}}{Matrix of true modal category values for 
#'   each exposure item and latent class. pxK}
#'   \item{\code{true_global_thetas}}{Array of true item level probabilities for 
#'   each exposure item and latent class. pxKxR}
#'   \item{\code{true_xi}}{Matrix of true probit regression model coefficients. Kxq}
#'   \item{\code{true_Ci}}{Vector of true latent class assignments. nx1}
#'   \item{\code{true_Phi_mat}}{Matrix of true outcome probabilities for each
#'   class and stratum. KxS}
#'   \item{\code{true_Phi}}{Vector of true outcome probabilities for all 
#'   individuals. nx1}
#' }
#' @source Simulated data using `simulate_data.R`
#' @keywords datasets
#' @examples 
#' data(sim_data)
#' x_mat <- sim_data$X_data
"sim_data"