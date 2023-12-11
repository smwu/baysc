#===================================================================
# Simulate data with stratum-specific class-membership probabilities
# Informative sampling: strata variable influences class and outcome
# Programmer: SM Wu
# Date updated: 2023/04/10
#===================================================================

# library(LaplacesDemon) # Sample from distributions
# library(SimCorMultRes) # Simulated correlated binary outcomes
# library(dplyr)         # Data wrangling and manipulation
# library(ggpubr)        # Multiple plots side-by-side


#' Get unique categorical probabilities for multinomial logistic regresion
#' 
#' @description
#' Get categorical probabilities for a unique design matrix corresponding to 
#' specified multinomial logistic regression formula and coefficient parameters
#' 
#' @param beta_mat Coefficient parameters for the linear predictor terms of a 
#' multinomial logistic regression. KxQ, where Q is the number of covariate 
#' terms in the regression
#' @param formula Formula for multinomial logistic regression. Should start with
#' `"~ c_all"` if generating exposure X. All variables must be found in `V`.
#' @param V Dataframe of variables 
#' @return
#' List `categ_probs` containing:
#' \describe{
#'   \item{\code{unique_comb}}{Matrix of unique covariate combinations in the
#'   design matrix created using the specified formula. mxQ, where m is the 
#'   number of unique covariate combinations.}
#'   \item{\code{unique_prob_mat}}{Matrix of category probabilities for all 
#'   unique covariate cominations. mxQ}
#' }
#' @importFrom stats model.matrix as.formula
#' @keywords internal
#' @export
#' @examples
#' # Get pi probabilities for C ~ S
#' V_unique <- data.frame(s_all = as.factor(1:2))
#' formula_lc <- "~ s_all"
#' beta_mat <- matrix(c(0, 0,  # Kx2 (first row is all 0's)
#'                      0.5, 1.3,
#'                      -0.4, 1.5), byrow = TRUE, nrow = 3, ncol = 2)
#' categ_probs_pi <- get_categ_probs(beta_mat = beta_mat, formula = formula_lc, 
#'                                V = V_unique)
#' 
#' # Get theta probabilities for X ~ C + S
#' V_unique <- expand.grid(factor(1:3), factor(1:2))  
#' colnames(V_unique) <- c("c_all", "s_all")    
#' formula <- "~ c_all + s_all"
#' non_div_mode <- log(0.05 / 0.85)
#' mode_div_non <- log(0.85 / 0.05)
#' # r=1: 0.85, 0.05, 0.05; None, High, Med
#' beta_mat <- matrix(c(0,            0,              0,              0,
#'                      non_div_mode, mode_div_non,   mode_div_non,   0,
#'                      non_div_mode, mode_div_non,   2*mode_div_non, 0,
#'                      non_div_mode, 2*mode_div_non, mode_div_non,   0),
#'                    nrow = 4, byrow = TRUE)
#' categ_probs_theta <- get_categ_probs(beta_mat = beta_mat, formula = formula, 
#'                                      V = V_unique)
#' 
#' # Add in association with S where categ High has prob 0.91 for S=1
#' beta_mat[-1, 4] <- c(0, 0, log(0.91/0.03) - mode_div_non)
#' categ_probs_theta_s <- get_categ_probs(beta_mat = beta_mat, formula = formula, 
#'                                        V = V_unique)
#' 
get_categ_probs <- function(beta_mat, formula, V) {
  # Create design matrix
  design_mat <- stats::model.matrix(stats::as.formula(formula), data = V)
  # Get unique covariate combinations in design matrix
  unique_comb <- unique(design_mat)
  
  # Get linear predictor terms for each unique combo 
  unique_lin_preds <- unique_comb %*% t(beta_mat)
  # Obtain multinomial category probabilities, normalized to sum to 1
  unique_numer <- exp(unique_lin_preds)
  unique_denom <- rowSums(unique_numer)
  unique_prob_mat <- unique_numer / unique_denom
  
  # Return unique design matrix covariate combos and multinomial category probs
  categ_probs <- list(unique_comb = unique_comb, 
                      unique_prob_mat = unique_prob_mat)
  return(categ_probs)
}

#' Obtain betas matrix for generating categorical latent class assignment 
#' 
#' @description
#' Obtain matrix of betas that can be used to generate the categorical latent 
#' class assignment variable C, depending on binary variable S, using the 
#' specified matrix of assignment probabilities, \eqn{\pi}, and multinomial 
#' logistic regression formula \eqn{C ~ S}
#' @param pi_mat Matrix of class assignment probabilities, pi. SxK
#' @param design_mat Design matrix containing s_all
#' @return 
#' Returns `beta_mat` matrix of betas to be used in a multinomial logistic 
#' regression to generate a categorical variable C. `beta_mat` has R rows and 
#' number of columns equal to the number of columns in the design matrix
#' @keywords internal
#' @export
#' @examples
#' pi_mat <- matrix(c(0.3, 0.5, 0.2, 
#'                    0.1, 0.6, 0.3), byrow = TRUE, nrow = 2)
#' formula_lc <- "~ s_all"
#' V <- data.frame(s_all = as.factor(c(rep(0, times = 60000),
#'                                     rep(1, times = 20000))))
#' design_mat <- model.matrix(as.formula(formula_lc), data = V)
#' get_betas_c(pi_mat = pi_mat, design_mat = design_mat)
get_betas_c <- function(pi_mat, design_mat) {
  K <- ncol(pi_mat)
  beta_mat <- matrix(0, nrow = K, ncol = ncol(design_mat))
  for (k in 2:K) {
    beta_mat[k, 1] <- log(pi_mat[1, k] / pi_mat[1, 1])
    beta_mat[k, 2] <- log(pi_mat[2, k] / pi_mat[2, 1]) - beta_mat[k, 1]
  }
  return(beta_mat)
}

#' Obtain beta matrices for generating multivariate categorical exposures
#' 
#' @description
#' Obtain list of beta matrices that can be used to generate the multivariate 
#' categorical exposure variable X using a multinomial logistic regression 
#' covariates categorical latent class C and, if desired, categorical variable S, 
#' using the specified matrix of exposure category probabilities, \eqn{\theta}, 
#' and an input design matrix. 
#' 
#' @param thetas Matrix of exposure category probabilities, \eqn{\theta}. JxK, 
#' where J is the number of exposure items and K is the number of latent classes.
#' @param modal_theta_prob Probability of true exposure level. Default is 0.85.
#' @param R Number of exposure levels. Fixed across exposure items.
#' @param design_mat Design matrix of unique combinations of covariate variable 
#' values, or full design matrix for all individuals.
#' @param depends_s Boolean specifying whether X should also have an association
#' with S. If true, level R (i.e., the last level) has an additional 0.02*(R-1)
#' probability added to it among those with S=H, where H is the number of levels 
#' of S. 
#' @return 
#' Returns list `beta_mat_x` of length J with each element a matrix of betas to 
#' be used in a multinomial logistic regression to generate a categorical 
#' exposure variable for that item. Each matrix of betas has R rows and number 
#' of columns equal to the number of columns in the design matrix.
#' @seealso [simulate_pop()]
#' @export
#' @examples
#' ## Example 1: X ~ C
#' # Function to get beta's from a set of theta's
#' J <- 30
#' thetas <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J), 
#'                                       rep(3, times = 0.5 * J)),
#'                                C2 = c(rep(4, times = 0.2 * J), 
#'                                       rep(2, times = 0.8 * J)),
#'                                C3 = c(rep(3, times = 0.3 * J), 
#'                                       rep(4, times = 0.4 * J),
#'                                       rep(1, times = 0.3 * J))))
#' modal_theta_prob <- 0.85
#' R <- 4
#' # Dataframe of unique values of covariates
#' V_unique <- data.frame(c_all = as.factor(1:3))
#' formula <- "~ c_all"  
#' # Design matrix of unique combinations of covariate variable values
#' design_mat_unique <- model.matrix(as.formula(formula), data = V_unique)
#' beta_mat_x <- get_betas_x(thetas = thetas, 
#'                           modal_theta_prob = modal_theta_prob, R = R, 
#'                           design_mat = design_mat_unique)
#' 
#' ## Example 2: X ~ C + S
#' J <- 30
#' thetas <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J), 
#'                                       rep(3, times = 0.5 * J)),
#'                                C2 = c(rep(4, times = 0.2 * J), 
#'                                       rep(2, times = 0.8 * J)),
#'                                C3 = c(rep(3, times = 0.3 * J), 
#'                                       rep(4, times = 0.4 * J),
#'                                       rep(1, times = 0.3 * J))))
#' # Dataframe of unique values of covariates                                   
#' V_unique <- expand.grid(factor(1:3), factor(1:2))  
#' colnames(V_unique) <- c("c_all", "s_all")                                
#' modal_theta_prob <- 0.85
#' R <- 4
#' formula <- "~ c_all + s_all"  
#' # Design matrix of unique combinations of covariate variable values
#' design_mat_unique <- model.matrix(as.formula(formula), data = V_unique)
#' beta_mat_x <- get_betas_x(thetas = thetas, 
#'                           modal_theta_prob = modal_theta_prob, R = R, 
#'                           design_mat = design_mat_unique, depends_s = TRUE)
#'                           
get_betas_x <- function(thetas, modal_theta_prob = 0.85, R, design_mat, 
                        depends_s = FALSE) {
  # Get dimensions and initialize values
  K <- ncol(thetas)
  J <- nrow(thetas)
  non_mode <- (1 - modal_theta_prob) / K
  mode_div_non <- log(modal_theta_prob / non_mode)
  non_div_mode <- log(non_mode / modal_theta_prob)
  beta_mat_x <- vector(mode = "list", length = J)
  Q <- ncol(design_mat)
  
  # For each exposure item j, get beta matrix
  for (j in 1:J) {
    beta_mat_j <- matrix(0, nrow = R, ncol = Q)
    theta_j <- thetas[j, ]  # thetas for item j
    # Case when mode for C=1 is category 1
    if (theta_j[1] == 1) {  
      # C = 1
      beta_mat_j[-1, 1] <- non_div_mode
      # C = 2, 3, ...
      for (k in 2:K) {
        # If mode is same as C=1 mode, then all 0's
        # Else if mode is different from C = 1
        if (theta_j[k] != theta_j[1]) { 
          # beta_rk = M + M*I(r is mode)
          beta_mat_j[-1, k] = mode_div_non
          beta_mat_j[theta_j[k], k] = beta_mat_j[theta_j[k], k] + mode_div_non
        }
      }
    # Case when mode for C=1 is not category 1
    } else {  
      # C = 1
      beta_mat_j[theta_j[1], 1] <- mode_div_non 
      # C = 2, 3, ...
      for (k in 2:K) {
        if (theta_j[k] != theta_j[1]) { # Mode is different from C = 1
          if (theta_j[k] == 1) {  # Mode is category 1
            # beta_rk = N + N*I(r is C=1 mode)
            beta_mat_j[-1, k] = non_div_mode
            beta_mat_j[theta_j[1], k] = beta_mat_j[theta_j[1], k] + non_div_mode
          } else {  # Mode is not category 1
            # beta_rk = N*I(r is C=1 mode) + M*I(r is mode)
            beta_mat_j[theta_j[1], k] = non_div_mode
            beta_mat_j[theta_j[k], k] = mode_div_non
          }
        }
      }
    }
    
    if (depends_s) {
      # If desired, add in association with S where categ R has prob 0.02*(R-1)
      # higher probability for those with S=H
      new_non_mode <- non_mode - 0.02
      new_mode <- 1 - (new_non_mode * (R - 1))
      beta_mat_j[-1, Q] <- c(0, 0, log(new_mode/new_non_mode) - mode_div_non)
    }
    
    # Add beta matrix to list of betas for all items
    beta_mat_x[[j]] <- beta_mat_j
  }
  # Return list of beta matrices
  return(beta_mat_x)
}



#' Create a categorical variable using multinomial logistic regression
#' 
#' @param beta_mat Coefficient parameters for the linear predictor terms of a 
#' multinomial logistic regression. KxQ, where K is the number of categories 
#' and Q is the number of covariate terms in the regression
#' @param design_mat Regression design matrix. NxQ
#' @param split_dim Optional covariate to get separate category probabilities 
#' for. Default is `NULL`. If not `NULL`, must also specify `V` containing the
#' covariate of interest.
#' @param V Optional dataframe of covariate variables. Default is `NULL`. If 
#' `split_dim` is not `NULL`, must be specified as a dataframe of dimension NxQ' 
#' where Q' is the number of covariate variables.
#' @return 
#' Returns list `out_vars` containing
#' \describe{
#'   \item{\code{categ_var}}{Generated categorical variable. Nx1}
#'   \item{\code{pi_mat}}{Matrix of category probabilities for all individuals. NxK}
#'   \item{\code{true_pi}}{Vector of overall category probabilities across all
#'   individuals. Nx1}
#'   \item{\code{pi_split}}{If `split_dim` is specified, HxK matrix of category 
#'   probabilities split by the levels in `split_dim`, where H is the number of 
#'   levels in the `split_dim` variable, and K is the number of categories as 
#'   specified by `beta_mat`. Otherwise, `NULL` if `split_dim` is `NULL`.}
#' }
#' @seealso [simulate_pop()]
#' @export
#' @examples
#' # Define multinomial logistic regression formula to depend on s_all
#' formula_c <- "~ s_all"
#' # Create dataframe of covariate variables
#' V <- data.frame(s_all = as.factor(c(rep(0, times = 60000), 
#'                                     rep(1, times = 20000))))
#' # Create design matrix
#' design_mat_c <- model.matrix(as.formula(formula_c), data = V) 
#' # Specify matrix of multinomial logistic regression coefficients for each
#' # category. Kx2 matrix with first row all 0's.
#' # The coefficients below corresponds to category probabilities of 
#' # (0.3, 0.5, 0.2) for S=1 and (0.1, 0.6, 0.3) for S=2, giving overall 
#' # probabilities (0.253, 0.522, 0.225)
#' beta_mat_c <- matrix(c(0, 0, 
#'                        0.5, 1.3,
#'                        -0.4, 1.5), byrow = TRUE, nrow = 3, ncol = 2)
#' # Create categorical variable                      
#' out_vars <- create_categ_var(beta_mat = beta_mat_c, design_mat = design_mat_c,
#'                              split_dim = "s_all", V = V)
#' out_vars$categ_var
create_categ_var <- function(beta_mat, design_mat, split_dim = NULL, V = NULL) {
  if (!is.matrix(design_mat)) {
    stop("design_mat must be a dataframe")
  }
  if (!is.matrix(beta_mat)) {
    stop("beta_mat must be a matrix")
  } else if (ncol(beta_mat) != ncol(design_mat)) {
    stop("beta_mat must have the same number of columns as design_mat")
  } else if (sum(beta_mat[1, ]) != 0) {
    stop("first row of beta_mat should be all 0's")
  }
  if (!is.null(split_dim)) {
    if (is.null(V)) {
      stop("V dataframe must be specified if split_dim is not NULL")
    } else if (!(split_dim %in% colnames(V))) {
      stop("split_dim must be in colnames(V)")
    } else if (!is.factor(V[, colnames(V) == split_dim])) {
      stop("split_dim must specify a factor variable in V")
    }
  }
  
  # Obtain dimensions and initialize matrices
  K <- nrow(beta_mat)   # Number of categories
  N <- nrow(design_mat) # Number of individuals
  lin_preds <- matrix(NA, nrow = N, ncol = K)
  pi_mat <- matrix(NA, nrow = N, ncol = K)
  
  # Obtain multinomial logistic regression linear predictors for each individual
  # and each category
  for (k in 1:K) {
    lin_preds[, k] <- design_mat %*% beta_mat[k, ] 
  }
  denom <- exp(lin_preds) %*% rep(1, K)
  # Obtain matrix of category probabilities for each individual
  for (k in 1:K) {
    pi_mat[, k] <- exp(lin_preds[, k]) / denom
  }
  # Normalize to sum to 1
  pi_mat <- pi_mat / rowSums(pi_mat)
        # temp <- apply(pi_mat, 1, sum)
        # sum(temp - 1)
        # all(temp == 1)
  
  # Create categorical variable
  categ_var <- sapply(1:N, function(i) sample(1:K, size = 1, prob = pi_mat[i, ]))
        # prop.table(table(categ_var))
        # prop.table(table(categ_var[1:60000]))
        # prop.table(table(categ_var[-(1:60000)]))
  
  # Overall category probabilities observed in the population
  true_pi <- prop.table(table(categ_var))
  
  # If matrix split on additional variable is desired, create and add to output
  if (!is.null(split_dim)) {
    split_var <- V[, colnames(V) == split_dim]
    n_levels <- length(levels(split_var))
    pi_split <- matrix(NA, nrow = n_levels, ncol = K)
    for (h in 1:n_levels) {
      pi_split[h, ] <- colMeans(pi_mat[split_var == h, ])
    }
  } else {
    pi_split <- NULL
  }
  
  # Return output
  out_vars <- list(categ_var = categ_var, pi_mat = pi_mat, true_pi = true_pi,
                   pi_split = pi_split)
  
  return(out_vars)
}

### Equivalent way for categorical variable
# pi_s <- matrix(c(0.3, 0.5, 0.2, 
#                  0.1, 0.6, 0.3), nrow = 2, byrow = TRUE)
# s_all <- V$s_all
# c_all <- numeric(N)
# for (i in 1:N) {
#   c_all[i] <- rcat(n = 1, p = pi_s[s_all[i], ])
# }


#' Create and save simulated population data
#' 
#' @description
#' `simulate_pop` creates and saves simulated population data according to input
#' specifications
#' 
#' @param N Population size. Default is 80000.
#' @param J Number of exposure items. Default is 30.
#' @param K Number of latent classes. Default is 3.
#' @param R Number of exposure categories for all items. Default is 4.
#' @param N_s Population size for each level of categorical variable S. 
#' Default is `c(60000, 20000)`, corresponding to a binary S with H=2 levels.
#' @param formula_c String specifying formula for multinomial logistic 
#' regression to create category latent class assignment C. Default is `"~ s_all"`, 
#' which generates C dependent on stratum variable S. All variables in the 
#' formula must be `"c_all"`, `"s_all"`, or specified in `V_additional`.
#' @param formula_x String specifying formula for multinomial logistic regression 
#' to create multivariate categorical exposure X. Default is `"~ c_all"`, 
#' which generates X dependent on latent class C. All variables in 
#' the formula must be `"c_all"`, `"s_all"`, or specified in `V_additional`.
#' @param formula_y String specifying formula for logistic regression to create 
#' binary outcome Y. Default is `"~ c_all * s_all"`, which generates Y dependent 
#' on latent class C, stratum S, and their interaction. All variables in the 
#' formula must be `"c_all"`, `"s_all"`, or specified in `V_additional`.
#' @param modal_theta_prob Probability between 0 and 1 for the most likely 
#' exposure category, assumed to be the same for all items. For all other
#' categories, 1 - `modal_theta_prob` is evenly split across them. Default is 
#' 0.85, so if there are R=4 exposure categories, then the 3 non-modal 
#' categories occur with probability \eqn{(1-0.85)/3 = 0.05}
#' @param beta_mat_c Coefficient parameters for a multinomial logistic 
#' regression to generate latent class assignment C. KxQ, where K is the 
#' number of categories and Q is the number of covariate terms in the regression. 
#' Default is `NULL` and default values are used. See details.
#' @param beta_mat_x Coefficient parameters for a multinomial logistic 
#' regression to generate multivariate exposure X. List of J matrices, 
#' each of dimension KxQ. Default is `NULL` and default values are used. See details.
#' @param beta_vec_y Coefficient parameters for a logistic regression to 
#' generate binary outcome Y. Qx1. Default is `NULL` and 
#' default values are used. See details.
#' @param xi_mat_y Alternative to `xi_mat_y` but in mixture reference coding. 
#' Default is `NULL`. If specified, must have K rows and q columns, where q is 
#' the number of covariate terms in the regression, excluding all terms 
#' involving C.
#' @param V_additional Dataframe with any additional variables other than C and
#' S to be used to generate variables. Number of rows should be N. Default is `NULL`.
#' @param cluster_size Size of clusters for clustering in the outcome, Y, assumed
#' to be the same for all clusters. Default is 80, corresponding to 1000 
#' clusters if N is 80000, assuming exchangeable latent correlation matrix with 
#' correlation 0.5 on the off-diagonals.
#' @param pop_seed Numeric seed for population generation. Default is 1.
#' @param save_res Boolean specifying if results should be saved. Default = `TRUE`.
#' @param save_path String specifying directory and file name to save results,
#' e.g., "~/Documents/stratified_classes". Default is `NULL`.
#' 
#' @details
#' If `NULL` (default) is used for `beta_mat_c`, `beta_mat_x`, or `beta_mat_y`, 
#' the following default values are used, corresponding to the scenario with K=3 
#' latent class and H=2 levels for variable S: `beta_mat_c` is a 3x2 matrix with
#' values `c(0, 0, 0.5, 1.3, -0.4, 1.5)` by row; `beta_mat_x` is a list of J=30
#' 4x4 matrices as in Example 1 provided below; and `beta_vec_y` is a vector 
#' of length \eqn{3*2=6} with values `c(1, -0.7, -1.5, -0.5, -0.5, -0.3)`. The 
#' generation of default values is demonstrated in Example 1 provided below.
#' 
#' To save the simulated population, set `save_res = TRUE` (default) and 
#' `save_path` to a string that specifies both the location and the beginning of 
#' the file name (e.g., "~/Documents/stratified_classes"). The file name will 
#' have "sim_pop.RData" appended to it.
#' 
#' @return
#' Returns list `sim_pop` containing:
#' \describe{
#'   \item{\code{N}}{Population size}
#'   \item{\code{J}}{Number of exposure items}
#'   \item{\code{R}}{Number of exposure categories}
#'   \item{\code{H}}{Number of stratum (i.e., levels of S)}
#'   \item{\code{N_s}}{Vector of population sizes for levels of S. Hx1}
#'   \item{\code{true_K}}{True number of latent classes, K}
#'   \item{\code{true_Ai}}{`NULL` or vector of additional continuous variable if 
#'   `a_all` is provided in `V_additional`. Nx1}
#'   \item{\code{true_Bi}}{`NULL` or vector of additional binary variable if 
#'   `b_all` is provided in `V_additional`. Nx1}
#'   \item{\code{true_Si}}{Vector of true stratum indicators. Nx1}
#'   \item{\code{true_Ci}}{Vector of true latent class indicators. Nx1}
#'   \item{\code{true_pi}}{Vector of true pi values overall in the population. Kx1}
#'   \item{\code{true_pi_s}}{`NULL` or HxK matrix of true pi values within each 
#'   level of S}
#'   \item{\code{X_data}}{Matrix of multivariate categorical exposure for all 
#'   individuals. NxJ}
#'   \item{\code{true_global_patterns}}{Matrix of true global exposure patterns 
#'   defined by modal category. JxK}
#'   \item{\code{true_global_thetas}}{Array of true thetas. JxKxR}
#'   \item{\code{Y_data}}{Vector of binary outcome for all individuals. Nx1}
#'   \item{\code{cluster_id}}{Vector of true cluster indicators. Nx1}
#'   \item{\code{cluster_size}}{Cluster size}
#'   \item{\code{true_xi}}{Matrix of probit regression coefficients in mixture 
#'   reference coding. Kxq, where q is the number of covariate terms in the 
#'   regression, excluding all terms involving C}
#'   \item{\code{true_Phi}}{Vector of true outcome probabilities for all 
#'   individuals. Nx1}
#'   \item{\code{true_Phi_mat}}{`NULL` or matrix of true outcome probabilities
#'   for individuals aggregated by C and S. KxH}
#' }
#' If `save_res = TRUE` (default), also saves `sim_pop` as 
#' `[save_path]_sim_pop.RData`. 
#' 
#' @seealso [simulate_samp()] [create_categ_var()] [get_betas_x()]
#' @importFrom stringr str_detect
#' @importFrom stats terms as.formula model.matrix rnorm rbinom pnorm toeplitz
#' @importFrom SimCorMultRes rbin
#' @export
#' @examples 
#' ### Example 1: Default values
#' # Default dimensions
#' N = 80000; J = 30; K = 3; R = 4; N_s = c(60000, 20000)
#' H <- length(N_s); modal_theta_prob = 0.85; cluster_size = 80; pop_seed <- 1
#' 
#' # Generate C ~ S
#' formula_c <- "~ s_all"
#' beta_mat_c <- matrix(c(0, 0, 0.5, 1.3, -0.4, 1.5), 
#'                      byrow = TRUE, nrow = K, ncol = H)
#' 
#' # Generate X ~ C
#' formula_x <- "~ c_all"
#' V_unique <- data.frame(c_all = factor(1:K)) # Df of unique cov values
#' design_mat_unique <- stats::model.matrix(stats::as.formula(formula_x), 
#'                                          data = V_unique)
#' thetas <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
#'                                       rep(3, times = 0.5 * J)),
#'                                C2 = c(rep(4, times = 0.2 * J),
#'                                       rep(2, times = 0.8 * J)),
#'                                C3 = c(rep(3, times = 0.3 * J),
#'                                       rep(4, times = 0.4 * J),
#'                                       rep(1, times = 0.3 * J))))
#' beta_mat_x <- get_betas_x(thetas = thetas, modal_theta_prob = modal_theta_prob,
#'                           R = R, design_mat = design_mat_unique)
#' 
#' # Generate Y ~ C + S + C:S
#' formula_y <- "~ c_all * s_all"
#' beta_vec_y <- c(1, -0.7, -1.5, -0.5, -0.5, -0.3)
#' 
#' # Create population
#' sim_pop <- simulate_pop(N = N, J = J, K = K, R = R, N_s = N_s,
#'                         modal_theta_prob = modal_theta_prob, 
#'                         formula_c = formula_c, formula_x = formula_x, 
#'                         formula_y = formula_y, beta_mat_c = beta_mat_c, 
#'                         beta_mat_x = beta_mat_x, beta_vec_y = beta_vec_y, 
#'                         cluster_size = cluster_size, 
#'                         pop_seed = pop_seed, save_res = FALSE)
#'                         
#' \dontrun{
#' ### Example 2: Selection-dependent pattern profiles     
#' # Generate X ~ C + S
#' formula_x <- "~ c_all + s_all"
#' V_unique <- expand.grid(factor(1:K), factor(1:H)) # Df of unique cov values
#' colnames(V_unique) <- c("c_all", "s_all")
#' design_mat_unique <- stats::model.matrix(stats::as.formula(formula_x), 
#'                                          data = V_unique)
#' beta_mat_x <- get_betas_x(thetas = thetas, modal_theta_prob = modal_theta_prob,
#'                           R = R, design_mat = design_mat_unique, depends_s = TRUE)
#' 
#' # Create population
#' sim_pop <- simulate_pop(N = N, J = J, K = K, R = R, N_s = N_s,
#'                         modal_theta_prob = modal_theta_prob, 
#'                         formula_c = formula_c, formula_x = formula_x, 
#'                         formula_y = formula_y, beta_mat_c = beta_mat_c, 
#'                         beta_mat_x = beta_mat_x, beta_vec_y = beta_vec_y, 
#'                         cluster_size = cluster_size, 
#'                         pop_seed = pop_seed, save_res = FALSE)
#' }
#' 
#' \dontrun{
#' ### Example 2: Additional effect modifiers for Y      
#' # Continuous variable A for age centered about 0
#' a_all <- stats::rnorm(n = N, mean = 0, sd = 5)
#' # Binary variable B for physically inactive or active
#' b_all <- as.factor(stats::rbinom(n = N, size = 1, prob = 0.3) + 1)
#' # Create dataframe of additional variables A and B
#' V_additional <- data.frame(a_all, b_all)    
#'               
#' # Generate Y ~ C + S + A + B + C:S + C:A + C:B
#' formula_y <- "~ c_all * (s_all + a_all + b_all)"
#' beta_vec_y <- c(1, -0.7, -1.5, -0.5, -0.5, -0.3, -0.04, 0.09, 0.08, 0.4, 
#'                 -0.7, -0.6)
#'                 
#' # Create population
#' sim_pop <- simulate_pop(N = N, J = J, K = K, R = R, N_s = N_s,
#'                         modal_theta_prob = modal_theta_prob, 
#'                         formula_c = formula_c, formula_x = formula_x, 
#'                         formula_y = formula_y, beta_mat_c = beta_mat_c, 
#'                         beta_mat_x = beta_mat_x, beta_vec_y = beta_vec_y, 
#'                         V_additional = V_additional, cluster_size = cluster_size, 
#'                         pop_seed = pop_seed, save_res = FALSE)
#' }
simulate_pop <- function(N = 80000, J = 30, K = 3, R = 4, 
                         N_s = c(60000, 20000),  modal_theta_prob = 0.85, 
                         formula_c = "~ s_all", 
                         formula_x = "~ c_all", 
                         formula_y = "~ c_all * s_all", 
                         beta_mat_c = NULL, beta_mat_x = NULL,
                         beta_vec_y = NULL, xi_mat_y = NULL, 
                         V_additional = NULL, 
                         cluster_size = 80, pop_seed = 1, 
                         save_res = TRUE, save_path = NULL) {

  # Set seed
  set.seed(pop_seed)
  #================== Check errors =============================================
  
  # Catch errors
  # Check modal_theta_prob is between 0 and 1
  if (modal_theta_prob < 0 | modal_theta_prob > 1) {
    stop("modal_theta_prob must be a probability between 0 and 1")
  }
  # Check V_additional is a dataframe
  if (!is.null(V_additional)) {
    if (!is.data.frame(V_additional)) {
      stop("V_additional must be a dataframe")
    }
  }
  # Check save_res and save_path
  if (!is.logical(save_res)) {
    stop("save_res must be a boolean specifying if results should be saved")
  } else if (save_res) {
    if (is.null(save_path) | !is.character(save_path)) {
      stop("save_path must be a string specifying a path and file name, such as 
           '~/Documents/stratified_classes'")
    } else {
      last_slash_ind <- regexpr("\\/[^\\/]*$", save_path)
      if (!dir.exists(substr(save_path, start = 1, stop = last_slash_ind))) {
        stop("directory specified in save_path does not exist")
      }
      if (last_slash_ind == length(save_path)) {
        stop("please append the start of a file name to the end of save_path. 
             For example, '~/Documents/stratified_classes' will produce a saved 
             file named 'stratified_classes_sim_pop.RData'")
      }
    }
  }
  # Check formulas for errors
  formulas <- c(formula_c, formula_x, formula_y)
  if (!all(is.character(formulas))) {
    stop("formula_c, formula_x, and formula_y must all be strings beginning 
    with '~' specifying formulas containing the variables that influence 
    latent class, exposure, and outcome, respectively")
  } else if (!(all(sapply(formulas, function(x) substring(x, 1, 1) == "~")))) {
    stop("formula_c, formula_x, and formula_y must all be strings beginning 
    with '~' specifying formulas containing the variables that influence 
    latent class, exposure, and outcome, respectively")
  }
  
  #================ Specify strength of dietary patterns =======================
  clust_mode <- modal_theta_prob
  non_mode <- (1 - clust_mode) / K
  
  #================ Create S variable ==========================================
  # Number of strata
  H <- length(N_s)
  # Create strata
  s_all <- unlist(sapply(1:H, function(x) rep(x, times = N_s[x])))
  # Dataframe of strata and additional variables
  if (!is.null(V_additional)) {
    V <- data.frame(s_all = as.factor(s_all), V_additional)
  } else {
    V <- data.frame(s_all = as.factor(s_all))
  }
  
  #================ Create latent class assignments C variable =================
  # Set defaults for generating latent class dependent on S and check beta_mat_c
  # for errors
  if (is.null(beta_mat_c)) {
    if (K != 3 | H != 2) {
      stop("If default for beta_mat_c is to be used, K must be equal to 3 and H 
           must be equal to 2.")
    }
    # Corresponds to pi = (0.3, 0.5, 0.2) for S=1 and (0.1, 0.6, 0.3) for S=2
    # for an overall true_pi ~= (0.253, 0.522, 0.225)
    beta_mat_c <- matrix(c(0, 0,  # Kx2 (first row is all 0's)
                         0.5, 1.3,
                         -0.4, 1.5), byrow = TRUE, nrow = K, ncol = H)
  } else if (!is.matrix(beta_mat_c)) {
    # Check that beta_mat_c is a matrix
    stop("beta_mat_c must be a matrix")
  }
  # Check that all variables in formula_c are available, removing interactions
  regr_vars <- labels(stats::terms(stats::as.formula(formula_c)))
  regr_vars <- regr_vars[!(stringr::str_detect(regr_vars, ":"))] 
  if (!all(regr_vars %in% colnames(V))) {
    stop("all variables in formula_c other than c_all and s_all must be 
      provided in V_additional")
  }
  
  # Create design matrix
  design_mat_c <- stats::model.matrix(stats::as.formula(formula_c), data = V)
  
  # Check congruency between beta_mat_c and formula_c
  if (ncol(beta_mat_c) != ncol(design_mat_c)) {
    stop(paste0("beta_mat_c must have columns corresponding to '", 
                paste0(colnames(design_mat_c), collapse = ", "),
                "' in the design matrix resulting from formula_c"))
  }
  
  # Create latent class assignments
  out_vars_c <- create_categ_var(beta_mat = beta_mat_c, 
                                 design_mat = design_mat_c, split_dim = "s_all",
                                 V = V)
  c_all <- out_vars_c$categ_var
  true_pi_s <- out_vars_c$pi_split
  true_pi <- out_vars_c$true_pi
  
  #================ Create multivariate exposure X variable ====================
  # Add latent class to data frame of variables
  V <- data.frame(c_all = as.factor(c_all), V)
  # Get design matrix
  design_mat_x <- stats::model.matrix(stats::as.formula(formula_x), data = V)
  
  # Set defaults for generating exposure X dependent on C and S and check
  # beta_mat_x for errors
  if (is.null(beta_mat_x)) {
    if (K != 3) {
      stop("If default for beta_mat_x is to be used, K must be equal to 3.")
    }
    thetas <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                          rep(3, times = 0.5 * J)),
                                   C2 = c(rep(4, times = 0.2 * J),
                                          rep(2, times = 0.8 * J)),
                                   C3 = c(rep(3, times = 0.3 * J),
                                          rep(4, times = 0.4 * J),
                                          rep(1, times = 0.3 * J))))
    beta_mat_x <- get_betas_x(thetas = thetas,
                              modal_theta_prob = modal_theta_prob, R = R,
                              design_mat = design_mat_x, depends_s = FALSE)
  } else if (!is.list(beta_mat_x) | !is.matrix(beta_mat_x[[1]])) {
    # Check that beta_mat_x is a list of matrices
    stop("beta_mat_x must be a list of matrices")
  } else if (ncol(beta_mat_x[[1]]) != ncol(design_mat_x)) {
    stop(paste0("Each matrix in beta_mat_x must have columns corresponding to '", 
                paste0(colnames(design_mat_x), collapse = ", "),
                "' in the design matrix resulting from formula_x"))
  }
  
  # Check that all variables in formula_x are available, removing interactions
  regr_vars <- labels(terms(as.formula(formula_x)))
  regr_vars <- regr_vars[!(stringr::str_detect(regr_vars, ":"))] 
  if (!all(regr_vars %in% colnames(V))) {
    stop("all variables in formula_x other than c_all and s_all must be 
    provided in V_additional")
  }
  
  # Create multivariate categorical exposure variable
  X_data <- matrix(NA, nrow = N, ncol = J)
  # Obtain underlying thetas for each item, class, and category
  true_global_thetas <- array(NA, dim = c(J, K, R))
  # Obtain underlying modal patterns for each item and class
  true_global_patterns <- matrix(NA, nrow = J, ncol = K)
  for (j in 1:J) {
    out_vars_x <- create_categ_var(beta_mat = beta_mat_x[[j]], 
                                   design_mat = design_mat_x, 
                                   split_dim = "c_all", V = V)
    X_data[, j] <- out_vars_x$categ_var
    true_global_thetas[j, , ] <- out_vars_x$pi_split
    true_global_patterns[j, ] <- apply(out_vars_x$pi_split, 1, which.max)
  }
  
  #================ Create binary Y outcome variable ===========================
  # Set defaults for generating exposure Y dependent on C and S and check
  # xi_mat_y and beta_vec_y for errors
  if (is.null(xi_mat_y)) {
    if (is.null(beta_vec_y)) {
      if (K != 3 | H != 2) {
        stop("If default for beta_vec_y is to be used, K must be equal to 3 and H 
           must be equal to 2.")
      }
      # Corresponds to 3x2 xi_mat_y with c(1, -0.5, 
      #                                    0.3, -1, 
      #                                    -0.5, -0.8) by row
      # and probabilities 0.84, 0.62, 0.31 for S=1 and 0.69, 0.24, 0.1 for S=2
      beta_vec_y <- c(1, -0.7, -1.5, -0.5, -0.5, -0.3)
      
    } else if (!is.vector(beta_vec_y)) {
      stop("beta_vec_y must be a vector if specified")
    } 
  } else {
    if (!is.matrix(xi_mat_y)) {
      # Check that xi_mat_y is a matrix
      stop("xi_mat_y must be a matrix")
    }
    
    # Convert regression coefficients from mixture reference to reference cell 
    beta_vec_y <- convert_mix_to_ref(est_xi = xi_mat_y)
  }
  
  # Check that all variables in formula_y are available, removing interactions
  regr_vars <- labels(stats::terms(stats::as.formula(formula_y)))
  regr_vars <- regr_vars[!(stringr::str_detect(regr_vars, ":"))] 
  if (!all(regr_vars %in% colnames(V))) {
    stop("all variables in formula_y other than c_all and s_all must be 
    provided in V_additional")
  }
  
  # Design matrix
  design_mat_y <- stats::model.matrix(stats::as.formula(formula_y), data = V)
  
  # Check design matrix and parameters congruence
  if (length(beta_vec_y) != ncol(design_mat_y)) {
    stop(paste0("Number of parameters specified in xi_mat_y or beta_vec_y must 
                match the number of parameters required for formula_y, which 
                produces a design matrix with columns '", 
                paste0(colnames(design_mat_y), collapse = ", "), "."))
  }
  
  # xi_mat_y is null, get corresponding values from beta_vec_y
  if (is.null(xi_mat_y)) {
    # xi is in mixture-reference coding (i.e., ref cell coding for each class/row):
    # xi1*I(C=1) + xi2*I(C=1,S=2) +
    # + xi3*I(C=2) + xi4*I(C=2,S=2)
    # + xi5*I(C=3) + xi6*I(C=3,S=2)
    # Number of covariate terms in the regression, excluding all terms with C
    cov_terms <- colnames(stats::model.matrix(stats::as.formula(formula_y)))
    cov_terms_no_c <- cov_terms[!stringr::str_detect(cov_terms, "c_all")]
    q <- length(cov_terms_no_c)
    xi_mat_y <- convert_ref_to_mix(K = K, q = q, est_beta = beta_vec_y)$est_xi
  }
  
  # Get vector of individual linear predictors
  lin_pred <- design_mat_y %*% beta_vec_y
  # Get vector of individual underlying outcome probabilities
  true_Phi_under <- stats::pnorm(lin_pred)
  
  ## Generate outcome Y depending on clustering
  # No clustering in the data
  if (cluster_size == 1) {
    cluster_id <- 1:N
    # Outcome data
    Y_data <- stats::rbinom(n = N, size = 1, prob = true_Phi_under)
    
  # Clustered data
  } else {
    # Modify formula_y to fit rbin() function specification
    formula_y_clus <- paste0("~ ", paste0(colnames(design_mat_y)[-1], 
                                          collapse = " + "))
    # Default is 1000 clusters of size 80
    # Simulate correlated binary outcomes 'SimCorMultRes' package
    # Assume exchangeable correlation matrix
    latent_correlation_matrix <- stats::toeplitz(c(1, rep(0.5, cluster_size - 1)))
    # Simulate correlated binary outcomes
    sim_binary <- SimCorMultRes::rbin(clsize = cluster_size, 
                                      intercepts = beta_vec_y[1],
                                      betas = beta_vec_y[-1], 
                                      xformula = formula_y_clus,
                                      xdata = design_mat_y,
                                      cor.matrix = latent_correlation_matrix, 
                                      link = "probit")
    # Cluster indicator for all individuals
    # By default, stratum 1 has clusters 1-250, stratum 2 has clusters 251-1000
    cluster_id <- sim_binary$simdata$id
    # Outcome data
    Y_data <- sim_binary$simdata$y
  }
  
  # If no additional confounders or interactions other than c_all*s_all, 
  # output outcome probabilities individually and aggregated by C and S
  y_covs <- labels(stats::terms(stats::as.formula(formula_y)))
  if (all(y_covs %in% c("c_all", "s_all", "c_all:s_all"))) {
    pop_inds <- 1:N
    # Matrix of finite population outcome probs aggregated by C and S
    true_Phi_mat <- matrix(NA, nrow=K, ncol=H)
    # Observed finite population probs for each individual based on C and S
    true_Phi <- numeric(N)
    for (h in 1:H) {
      for (k in 1:K) {
        # Individuals in stratum s class k
        N_s_k <- pop_inds[s_all == h & c_all == k]
        true_Phi_mat[k, h] <- sum(Y_data[N_s_k]) / length(N_s_k)
        true_Phi[N_s_k] <- true_Phi_mat[k, h]
      }
    }
    # If additional confounders and interactions, only output underlying 
    # individual outcome probabilities
  } else {
    true_Phi <- true_Phi_under
    true_Phi_mat <- NULL
  }
  
  #================ Save and return data =======================================
  # Define additional covariates for output
  true_Ai <- true_Bi <- NULL
  if ("a_all" %in% colnames(V)) {
    true_Ai <- V$a_all
  }
  if ("b_all" %in% colnames(V)) {
    true_Bi <- V$b_all
  }
  
  sim_pop <- list(N = N, J = J, R = R, H = H, N_s = N_s, true_K = K,
                  true_Ai = true_Ai, true_Bi = true_Bi, 
                  true_Si = s_all, true_Ci = c_all, true_pi = true_pi,
                  true_pi_s = true_pi_s, X_data = X_data, 
                  true_global_patterns = true_global_patterns,
                  true_global_thetas = true_global_thetas, Y_data = Y_data,
                  cluster_id = cluster_id, cluster_size = cluster_size,
                  true_xi = xi_mat_y, true_Phi = true_Phi, 
                  true_Phi_mat = true_Phi_mat)
  if (save_res) {
    save(sim_pop, file = paste0(save_path, "sim_pop.RData"))
  }

  return(sim_pop)
}



#' Create and save simulated sample data 
#' 
#' @description
#' `simulate_samp` creates and saves simulated sample data by sampling from the 
#' simulated population according to input specifications.
#' 
#' @inheritParams simulate_pop
#' @param sim_pop Simulated population output from `sim_pop()` function
#' @param samp_prop Proportion of population to sample. Default is 0.05. 
#' Either `samp_prop` or `samp_size` must be specified.
#' @param samp_size Sample size. Default is `NULL`. Either `samp_prop` or 
#' `samp_size` must be specified.
#' @param strat Boolean specifying whether to perform stratification by S.
#' Default is `TRUE`.
#' @param strat_dist Vector of relative stratum sizes in the sample. Components 
#' must be proportions that sum to 1. Default is `c(0.5, 0.5)`, which results in
#' equal stratum sizes in the sample. 
#' @param clust Boolean specifying whether to perform cluster sampling. Default 
#' is `TRUE`.
#' @param samp_seed Numeric seed for population generation. Default is 1.
#' 
#' @return
#' Returns list `sim_samp` containing:
#' \describe{
#'   \item{\code{samp_ind}}{Vector of population indices for sampled 
#'   individuals. nx1}
#'   \item{\code{sample_wt}}{Vector of sampling weights for sampled 
#'   individuals. nx1}
#'   \item{\code{N}}{Population size}
#'   \item{\code{J}}{Number of exposure items}
#'   \item{\code{R}}{Number of exposure categories}
#'   \item{\code{H}}{Number of stratum (i.e., levels of S)}
#'   \item{\code{N_s}}{Vector of population sizes for levels of S. Hx1}
#'   \item{\code{true_K}}{True number of latent classes, K}
#'   \item{\code{true_Ai}}{`NULL` or vector of additional continuous variable 
#'   for sampled individuals. nx1}
#'   \item{\code{true_Bi}}{`NULL` or vector of additional binary variable for 
#'   sampled individuals. nx1}
#'   \item{\code{true_Si}}{Vector of true stratum indicators for sampled 
#'   individuals. nx1}
#'   \item{\code{true_Ci}}{Vector of true latent class indicators for sampled 
#'   individuals. nx1}
#'   \item{\code{true_pi}}{Vector of true pi values overall in the population. Kx1}
#'   \item{\code{true_pi_s}}{`NULL` or HxK matrix of true pi values within each 
#'   level of S}
#'   \item{\code{X_data}}{Matrix of multivariate categorical exposure for 
#'   sampled individuals. nxJ}
#'   \item{\code{true_global_patterns}}{Matrix of true global exposure patterns 
#'   defined by modal category. JxK}
#'   \item{\code{true_global_thetas}}{Array of true thetas. JxKxR}
#'   \item{\code{Y_data}}{Vector of binary outcome for sampled individuals. nx1}
#'   \item{\code{cluster_id}}{Vector of true cluster indicators for sampled 
#'   individuals. nx1}
#'   \item{\code{cluster_size}}{Cluster size}
#'   \item{\code{true_xi}}{Matrix of probit regression coefficients in mixture 
#'   reference coding. Kxq, where q is the number of covariate terms in the 
#'   regression, excluding all terms involving C}
#'   \item{\code{true_Phi}}{Vector of true outcome probabilities for sampled 
#'   individuals. nx1}
#'   \item{\code{true_Phi_mat}}{`NULL` or matrix of true outcome probabilities
#'   for individuals aggregated by C and S. KxH}
#' }
#' If `save_res = TRUE` (default), also saves `sim_samp` as 
#' `[save_path]_sim_samp.RData`. 
#' 
#' @seealso [simulate_pop()] 
#' @export
#' @examples 
#' ### Example 1: Default values
#' # Create population
#' sim_pop <- simulate_pop(save_res = FALSE)
#' sim_samp <- simulate_samp(sim_pop = sim_pop, samp_prop = 0.05, strat = TRUE,
#'                           strat_dist = c(0.5, 0.5), clust = TRUE,
#'                           samp_seed = 101, save_res = FALSE)
simulate_samp <- function(sim_pop, samp_prop = 0.05, samp_size = NULL, 
                          strat = TRUE, strat_dist = c(0.5, 0.5), clust = TRUE, 
                          samp_seed = 101, save_res = TRUE, save_path) {
  # Set seed
  set.seed(samp_seed)
  #================ Check errors ===============================================
  # Check clustering
  if (sim_pop$cluster_size == 1) {
    warning("cluster sampling is requested but population does not contain 
                clustered data")
  }
  # Check stratification distribution
  if (strat) {
    if (is.null(strat_dist)) {
      stop("strat_dist must be specified if strat is TRUE")
    } else if (sum(strat_dist) != 1) {
      stop("strat_dist must be a vector of length H whose elements sum to 1")
    }
  }
  # Check save_res and save_path
  if (!is.logical(save_res)) {
    stop("save_res must be a boolean specifying if results should be saved")
  } else if (save_res) {
    if (is.null(save_path) | !is.character(save_path)) {
      stop("save_path must be a string specifying a path, such as '~/Documents/'")
    } else if (!dir.exists(save_path)) {
      stop("directory specified in save_path does not exist")
    }
  }
  
  #================ Get sample size ============================================
  if (is.null(samp_size)) {
    if (is.null(samp_prop)) {
      stop("one of samp_prop and samp_size must be specified")
    } else if (!(samp_prop > 0 & samp_prop <= 1)) {
      stop("samp_prop must be greater than 0 and less than or equal to 1")
    } else {
      samp_size <- ceiling(samp_prop * sim_pop$N)
    }
  } else if (!(samp_size > 0 & samp_size <= sim_pop$N)) {
    stop("samp_size must be greater than 0 and less than or equal to population 
         size N")
  }
  
  #================ Specify sampling design ====================================
  
  pop_inds <- 1:sim_pop$N  # Population indices
  cluster_id <- sim_pop$cluster_id  # Population cluster indicators
  sample_wt_temp <- numeric(sim_pop$N)  # Temp weights for population
  samp_ind_temp <- numeric(sim_pop$N)  # Temp sample indicator for population
  
  # Survey design is SRS
  if (!strat & !clust) {
    # Sampled individuals
    samp_ind_SRS <- sample(pop_inds, samp_size)
    samp_ind_temp[samp_ind_SRS] <- 1
    sample_wt_temp[samp_ind_SRS] <- sim_pop$N / samp_size
    
  # Survey design is stratified sampling with unequal probabilities
  } else if (strat & !clust) {
    if (length(strat_dist) != sim_pop$H) {
      stop(paste0("strat_dist must be a vector of length H = ", sim_pop$H, 
                  ", to match the number of levels in stratifying variable S"))
    }
    # Get stratum sample sizes
    n_s <- samp_size * strat_dist
    for (h in 1:sim_pop$H) {
      # Individuals in stratum h
      pop_s <- pop_inds[sim_pop$true_Si == h]
      # Sample from stratum h
      samp_ind_s <- sample(pop_s, n_s[h])
      samp_ind_temp[samp_ind_s] <- 1
      sample_wt_temp[samp_ind_s] <- sim_pop$N_s[h] / n_s[h]
    }
    
  # Survey design is one-stage cluster sampling
  } else if (clust & !strat) {
      # Number of clusters to sample
      n_clus <- samp_size / sim_pop$cluster_size
      # Check number of clusters is a whole number
      if (n_clus %% 1 == 0) {
        stop(paste0("for cluster sampling, samp_size must be a multiple of the cluster 
             size: ", sim_pop$cluster_size))
      }
      # Clusters
      clusters <- unique(sim_pop$cluster_id)
      # Sample clusters
      samp_clus_ind <- sample(clusters, size = n_clus)
      # Sampled individuals
      samp_ind <- pop_inds[sim_pop$cluster_id %in% samp_clus_ind]
      samp_ind_temp[samp_ind] <- 1
      sample_wt_temp[samp_ind] <- sim_pop$N / samp_size
    
  # Survey design is stratified cluster sampling
  } else {
    # Get stratum sample sizes
    n_s <- samp_size * strat_dist
    # Number of clusters to sample per stratum
    n_clus_s <- n_s / sim_pop$cluster_size
    # Check number of clusters are whole numbers
    if (!all(n_clus_s %% 1 == 0)) {
      stop(paste0("for stratified cluster sampling, samp_size and strat_dist 
                must be specified such that the sample size for each stratum 
                is a multiple of the cluster size: ", sim_pop$cluster_size, 
                  ". Currently, the strata are of size ", paste0(n_s, collapse = ", "), ". "))
    }
    for (h in 1:sim_pop$H) {
      # Clusters in stratum h 
      clus_s <- unique(sim_pop$cluster_id[sim_pop$true_Si == h])
      # Sample clusters from stratum h
      samp_clus_s <- sample(clus_s, size = n_clus_s[h])
      # Sampled individuals
      samp_ind_s <- pop_inds[sim_pop$cluster_id %in% samp_clus_s]
      # Update temp weights and sample indicators for population
      samp_ind_temp[samp_ind_s] <- 1
      sample_wt_temp[samp_ind_s] <- sim_pop$N_s[h] / n_s[h]
    }
  } 
  
  #================ Subset to sampled individuals ==============================
  samp_ind <- pop_inds[samp_ind_temp > 0]
  sample_wt <- sample_wt_temp[samp_ind]
  X_data <- sim_pop$X_data[samp_ind, ]
  Y_data <- sim_pop$Y_data[samp_ind]
  cluster_id <- sim_pop$cluster_id[samp_ind]
  true_Si <- sim_pop$true_Si[samp_ind]
  true_Ci <- sim_pop$true_Ci[samp_ind]
  true_Ai <- sim_pop$true_Ai[samp_ind]
  true_Bi <- sim_pop$true_Bi[samp_ind]
  true_Phi <- sim_pop$true_Phi[samp_ind]
  
  #================ Save and return output =====================================
  sim_samp <- list(samp_ind = samp_ind, sample_wt = sample_wt, 
                   N = sim_pop$N, J = sim_pop$J, R = sim_pop$R, H = sim_pop$H, 
                   N_s = sim_pop$N_s, true_K = sim_pop$true_K,
                   true_Ai = true_Ai, true_Bi = true_Bi, true_Si = true_Si,
                   true_Ci = true_Ci,  true_pi = sim_pop$true_pi, 
                   true_pi_s = sim_pop$true_pi_s, X_data = X_data, 
                   true_global_patterns = sim_pop$true_global_patterns,
                   true_global_thetas = sim_pop$true_global_thetas,
                   Y_data = Y_data, cluster_id = cluster_id, 
                   cluster_size = sim_pop$cluster_size,
                   true_xi = sim_pop$true_xi, true_Phi = true_Phi,
                   true_Phi_mat = sim_pop$true_Phi_mat)
  if (save_res) {
    save(sim_samp, file = paste0(save_path, "sim_samp.RData"))
  }
  return(sim_samp)
}




#==================== Additional Sanity checks =================================
# table(true_Si)
# prop.table(table(true_Ci))
# all(apply(true_global_thetas, c(1, 2), sum) == 1)
# apply(X_data[true_Ci == 1, ], 2, median)
# apply(X_data[true_Ci == 1, ], 2, mean)
# for (k in 1:K) {print(mean(Y_data[true_Ci == k]))}
# for (s in 1:S) {print(mean(Y_data[true_Si == s]))}
# sum(beta[c(1,2)]) == xi_vec[2]
# sum(beta[c(1,2,4,5)]) == xi_vec[5]
# sum(beta[c(1,3,4,6)]) == xi_vec[6]
# mean(sim_binary$simdata[sim_binary$simdata$C2S2 == 1, "y"])
# mean(sim_binary$simdata[sim_binary$simdata$C3S2 == 1, "y"])
# mean(sim_binary$simdata[sim_binary$simdata$S2 == 1, "y"])
# table(sim_binary$simdata$time)

# # Check that no cluster is in both strata
# intersect(unique(cluster_id[true_Si == 1]), unique(cluster_id[true_Si == 2]))
# table(sim_data$sample_wt)
# prop.table(table(sim_data$true_Si))
# prop.table(table(sim_data$true_Ci))

# # Check true outcome probabilities
# prop.table(table(round(true_Phi, 2)))
# mean(Y_data)

# # Check correlations between variables
# chisq.test(table(true_Ai, true_Ci))
# chisq.test(table(true_Ai, true_Si))
# chisq.test(table(true_Si, true_Ci))
# cor(true_Bi, true_Ci)

# # Plot density of B
# x_temp <- seq(-10,10,by=0.01)
# hist(true_Bi, breaks=30, freq = FALSE)
# lines(x_temp, dnorm(x_temp, mean=-5, sd=2)*(N_s[1]/N), col="red")
# lines(x_temp, dnorm(x_temp, mean=5, sd=2)*(N_s[2]/N), col="blue")

#================== Miscellaneous Old Code =====================================
### Code used to test `get_betas()` and `get_categ_probs()`
# # r=1: 0.05, 0.05, 0.05; Med, Low, Med
# beta_mat <- matrix(c(0,            0,              0,              0,  
#                      0,            mode_div_non,   0,              0,
#                      mode_div_non, non_div_mode,   0,              0,
#                      0,            0,              0,              0), 
#                    nrow = 4, byrow = TRUE)
# # r=1: 0.85, 0.85, 0.05; None, None, Med
# beta_mat <- matrix(c(0,            0,              0,              0,  
#                      non_div_mode, 0,              mode_div_non,   0,
#                      non_div_mode, 0,              2*mode_div_non, 0,
#                      non_div_mode, 0,              mode_div_non,   0), 
#                    nrow = 4, byrow = TRUE)
# # r=1: 0.05, 0.05, 0.85; Low, Low, None
# beta_mat <- matrix(c(0,            0,              0,              0,  
#                      mode_div_non, 0,              2*non_div_mode, 0,
#                      0,            0,              non_div_mode,   0,
#                      0,            0,              non_div_mode,   0), 
#                    nrow = 4, byrow = TRUE)
#
# # Population sanity checks
# prop.table(table(sim_pop$true_Ci[sim_pop$true_Si == 1]))
# prop.table(table(sim_pop$true_Ci[sim_pop$true_Si == 2]))
# prop.table(table(sim_pop$true_Ci))
# mean(sim_pop$Y_data)
# sim_pop$true_Phi_mat
# hist(sim_pop$true_Phi, breaks = 30)
# sim_pop$true_global_patterns
# 
# # Sample sanity Checks
# prop.table(table(sim_samp$true_Si))
# prop.table(table(sim_samp$true_Ci))
# prop.table(table(sim_samp$true_Ci[sim_samp$true_Si == 1]))
# prop.table(table(sim_samp$true_Ci[sim_samp$true_Si == 2]))
# K <- 3
# S <- 2
# samp_Phi_mat <- matrix(NA, nrow=K, ncol=S)
# for (k in 1:K) {
#   for (s in 1:S) {
#     samp_Phi_mat[k, s] <- sum(sim_samp$Y_data==1 & sim_samp$true_Si==s &
#                                 sim_samp$true_Ci==k) /
#       sum(sim_samp$true_Si==s & sim_samp$true_Ci==k)
#   }
# }
# samp_Phi_mat
# prop.table(table(sim_samp$X_data[sim_samp$true_Ci == 3,1] == 4))
# length(sim_samp$true_Phi)
# 
# #==================== Plot simulated theta modes ===============================
# # plot true modal theta patterns
# # Input:
# #   est_item_probs: JxKxR numeric matrix of true theta probabilities
# #   x_lab: String specifying x-axis label
# # Output: Plot of modal theta values for the K latent classes
# plot_sim_theta <- function(est_item_probs, x_lab) {
#   mode_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), which.max))
#   food_items <- 1:30
#   class_names <- 1:(dim(est_item_probs)[2])
#   rownames(mode_item_probs) <- food_items
#   colnames(mode_item_probs) <- class_names
#   mode_item_probs$Item <- rownames(mode_item_probs)
#   mode_plot <- mode_item_probs %>% gather("Class", "Level", -Item)
#   mode_plot %>% ggplot(aes(x=Class, y=factor(Item, levels = rev(food_items)),
#                            fill=factor(Level))) +
#     geom_tile(color="black", linewidth = 0.3) +
#     scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
#                       name = "Consumption Level",
#                       labels = c("None", "Low", "Med", "High")) +
#     xlab(x_lab) + ylab("Item") +
#     theme_classic() +
#     theme(text = element_text(size = 15),
#           axis.text.x = element_text(size = 11, color = "black"),
#           axis.text.y = element_text(size = 11, color = "black"),
#           axis.title.x = element_text(size = 13, color = "black", face = "bold"),
#           axis.title.y = element_text(size = 13, color = "black", face = "bold"),
#           legend.title = element_text(size = 13, color = "black", face = "bold"),
#           legend.text = element_text(size = 11, color = "black"),
#           legend.position = "right", legend.box.spacing = unit(0.2, "pt"))
# }
# 
# # Default setting
# samp_data_path <- paste0(wd, data_dir, "simdata_scen", 111211, "_iter",
#                          1, "_samp", 1, ".RData")
# load(samp_data_path)
# est_item_probs_default <- sim_data$true_global_thetas
# p_default <- plot_sim_theta(est_item_probs = est_item_probs_default,
#                             x_lab = "Default Pattern")
# 
# # Supervised setting
# samp_data_path <- paste0(wd, data_dir, "simdata_scen", 121211, "_iter",
#                          1, "_samp", 1, ".RData")
# load(samp_data_path)
# est_item_probs_supervised <- sim_data$true_global_thetas
# p_supervised <- plot_sim_theta(est_item_probs = est_item_probs_supervised,
#                                x_lab = "Overlap Pattern")
# 
# # Plot thetas together
# ggarrange(p_default, p_supervised, ncol = 2, common.legend = TRUE)
# p_default + theme(legend.position = "none",
#                   text = element_text(size = 15),
#                   axis.text.x = element_text(size = 13, color = "black"),
#                   axis.text.y = element_text(size = 13, color = "black"),
#                   axis.title.x = element_text(size = 15, color = "black", face = "bold"),
#                   axis.title.y = element_text(size = 15, color = "black", face = "bold"))
# p_supervised + theme(legend.position = "right",
#                      text = element_text(size = 15),
#                      axis.text.x = element_text(size = 13, color = "black"),
#                      axis.text.y = element_text(size = 13, color = "black"),
#                      axis.title.x = element_text(size = 15, color = "black", face = "bold"),
#                      axis.title.y = element_text(size = 15, color = "black", face = "bold"),
#                      legend.title = element_text(size = 15, color = "black", face = "bold"),
#                      legend.text = element_text(size = 15, color = "black"))
# 
