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
#' @param beta_mat Matrix of coefficient parameters for the linear predictor 
#' terms of a multinomial logistic regression for generating a categorical 
#' variable. Number of rows is equal to the number of levels in the categorical
#' variable. Number of columns is equal to the number of covariate terms in the 
#' regression.
#' @param formula Formula for multinomial logistic regression. Should start with
#' `"~ c_all"` if generating exposure X_i. All variables must be found in `V_unique`.
#' @param V_unique Dataframe with containing the unique values of the variables 
#' specified in `formula`.
#' @return
#' Matrix `categ_probs` where the first columns are the values of the covariates 
#' for the multinomial logistic regression, and the remaining columns are the 
#' corresponding probability of each category of the categorical outcome variable.  
#' 
#' @importFrom stats model.matrix as.formula
#' @export
#' @examples
#' # Get pi probabilities for C_i ~ S_i
#' K <- 3; S <- 2
#' formula_c <- "~ s_all"
#' V_unique <- data.frame(s_all = as.factor(1:S))
#' pi_mat <- matrix(c(0.3, 0.5, 0.2,   # class membership probs for S_i=1
#'                    0.1, 0.6, 0.3),  # class membership probs for S_i=2
#'                  byrow = TRUE, nrow = S, ncol = K)
#' beta_mat_c <- get_betas_c(pi_mat = pi_mat, formula_c = formula_c, 
#'                           V_unique = V_unique)
#' categ_probs_pi <- get_categ_probs(beta_mat = beta_mat_c, formula = formula_c, 
#'                                   V_unique = V_unique)
#' categ_probs_pi
#' 
#' # Get theta probabilities for X_i ~ C_i 
#' J <- 30; R <- 4; K <- 3
#' formula_x <- "~ c_all"
#' V_unique <- data.frame(c_all = as.factor(1:K))
#' profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J), 
#'                                         rep(3, times = 0.5 * J)),
#'                                  C2 = c(rep(4, times = 0.2 * J), 
#'                                         rep(2, times = 0.8 * J)),
#'                                  C3 = c(rep(3, times = 0.3 * J), 
#'                                         rep(4, times = 0.4 * J),
#'                                         rep(1, times = 0.3 * J))))
#' modal_theta_prob <- 0.85
#' beta_list_x <- get_betas_x(profiles = profiles, R = R, 
#'                            modal_theta_prob = modal_theta_prob, 
#'                            formula_x = formula_x, V_unique = V_unique)
#' categ_probs_theta <- get_categ_probs(beta_mat = beta_list_x[[1]], 
#'                                      formula = formula_x, V_unique = V_unique)
#' categ_probs_theta                                     
#' 
#' # Get theta probabilities for X_i ~ C_i + S_i
#' formula_x <- "~ c_all + s_all"
#' beta_list_x <- lapply(1:J, function(j) cbind(beta_list_x[[j]], 
#'                                              s_all = c(0, 0.5, 0, 0)))
#' V_unique <- expand.grid(c_all = as.factor(1:K), s_all = as.factor(1:S))                                            
#' categ_probs_theta_s <- get_categ_probs(beta_mat = beta_list_x[[1]], 
#'                                        formula = formula_x, V_unique = V_unique)
#' categ_probs_theta
#' 
get_categ_probs <- function(beta_mat, formula, V_unique) {
  # Create design matrix
  design_mat <- stats::model.matrix(stats::as.formula(formula), data = V_unique)
  # Get unique covariate combinations in design matrix
  unique_comb <- unique(design_mat)
  
  # Get linear predictor terms for each unique combo 
  unique_lin_preds <- unique_comb %*% t(beta_mat)
  # Obtain multinomial category probabilities, normalized to sum to 1
  unique_numer <- exp(unique_lin_preds)
  unique_denom <- rowSums(unique_numer)
  unique_prob_mat <- unique_numer / unique_denom
  if (!is.null(colnames(unique_prob_mat))) {
    colnames(unique_prob_mat) <- paste0("prob ", colnames(unique_prob_mat))
  }
  
  # Return multinomial category probs for the unique design matrix covariates
  categ_probs <- cbind(unique_comb[, -1], unique_prob_mat)

  return(categ_probs)
}

#' Obtain betas matrix for generating categorical latent class assignment 
#' 
#' @description
#' Obtain matrix of betas that can be used to generate the categorical latent 
#' class assignment variable C_i using a multinomial logistic regression where C_i 
#' may depend on a categorical covariate such as the stratum variable S_i. 
#' 
#' @param pi_mat Matrix where each row is the class membership probabilities for 
#' a level of the categorical covariate. SxK, where S is the number of levels of 
#' the categorical covariate and K is the number of latent classes. Rows of 
#' `pi_mat` must sum to 1.
#' @param formula_c String specifying formula for multinomial logistic 
#' regression to create category latent class assignment C_i.
#' @param V_unique Dataframe with one column containing the unique values of 
#' the categorical covariate specified in `formula_c`. If `formula_c = "~1"`, 
#' set `V_unique = NULL`.
#' 
#' @return 
#' Returns `beta_mat` matrix of betas to be used in a multinomial logistic 
#' regression to generate a categorical variable C_i. `beta_mat` has K rows and 
#' number of columns equal to the number of levels in the categorical covariate.
#' 
#' @importFrom stats terms as.formula model.matrix
#' @export
#' @examples
#' ## Example 1: latent class C_i depends on stratum variable S_i
#' # Number of latent classes and number of levels of S_i
#' K <- 3; S <- 2
#' # Formula specifying that C_i depends on S_i
#' formula_c <- "~ s_all"
#' # Dataframe with unique values of S_i
#' V_unique <- data.frame(s_all = as.factor(1:S))
#' # Matrix of class membership probabilities for each level of S_i
#' pi_mat <- matrix(c(0.3, 0.5, 0.2,   # class membership probs for S_i=1
#'                    0.1, 0.6, 0.3),  # class membership probs for S_i=2
#'                  byrow = TRUE, nrow = S, ncol = K)
#' # Get matrix of betas for generating C_i
#' beta_mat_c <- get_betas_c(pi_mat = pi_mat, formula_c = formula_c, 
#'                           V_unique = V_unique)
#' beta_mat_c
#' 
#' ## Example 2: latent class is generated independently of other variables
#' # Matrix of class membership probabilities
#' pi_mat <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
#' formula_c <- "~ 1"
#' V_unique <- NULL
#' get_betas_c(pi_mat = pi_mat, formula_c = formula_c, V_unique = V_unique)
#' 
get_betas_c <- function(pi_mat, formula_c, V_unique) {
  # Check errors
  if (!all(rowSums(pi_mat) == 1)) {
    stop("rows of pi_mat must sum to 1")
  }
  var_terms <- labels(stats::terms(stats::as.formula(formula_c)))
  # Only one variable allowed
  if (length(var_terms) > 1) {
    stop("formula_c cannot have more than one variable")
  }
  if (any(grepl(":", var_terms))) {
    stop("please do not specify interactions in formula_c")
  }
  
  # Number of latent classes
  K <- ncol(pi_mat)
  
  # If C_i is independent of other variables
  if (length(var_terms) == 0) {
    # Get beta matrix
    beta_mat <- matrix(0, nrow = K, ncol = 1)
    for (k in 2:K) {
      beta_mat[k, 1] <- log(pi_mat[1, k] / pi_mat[1, 1])
    }
    colnames(beta_mat) <- "(Intercept)"
    rownames(beta_mat) <- 1:K
  # If C_i depends on other variables
  } else {
    # Variable must be categorical and found in V_unique
    if (!is.factor(V_unique[[var_terms]])) {
      stop("V_unique must contain the unique values of a factor variable in formula_c")
    }
    # Create design matrix with unique values
    design_mat_unique <- stats::model.matrix(stats::as.formula(formula_c), 
                                             data = V_unique)
    # Number of levels of categorical covariate 
    S <- ncol(design_mat_unique)
    # Get beta matrix
    beta_mat <- matrix(0, nrow = K, ncol = S)
    for (k in 2:K) {
      beta_mat[k, 1] <- log(pi_mat[1, k] / pi_mat[1, 1])
      for (h in 2:S) {
        beta_mat[k, h] <- log(pi_mat[h, k] / pi_mat[h, 1]) - beta_mat[k, 1]
      }
    }
    colnames(beta_mat) <- colnames(design_mat_unique)
    rownames(beta_mat) <- 1:K
  }
  
  return(beta_mat)
}

#' Obtain list of beta matrices for generating multivariate categorical exposure X_i
#' 
#' @description
#' Obtain list of beta matrices that can be used to generate the multivariate 
#' categorical exposure variable X_i using a multinomial logistic regression where 
#' X_i depends on categorical covariate composed of latent class assignment C_i. 
#' 
#' @param profiles Matrix where each column is a latent class pattern profile 
#' and each row is the item level for all classes. JxK, where J is the number of 
#' exposure items and K is the number of latent classes. The item levels must 
#' range from 1 to R, where R is the number of levels for all items. 
#' @param R Number of exposure levels. Fixed across exposure items.
#' @param modal_theta_prob Probability of true exposure level. Default is 0.85.
#' @param formula_x String specifying formula for multinomial logistic 
#' regression to create multivariate categorical exposure X_i.
#' @param V_unique Dataframe with one column containing the unique values of 
#' the categorical covariate specified in `formula_x`. 
#' 
#' @return 
#' Returns list `beta_list_x` of length J with each element a matrix of betas to 
#' be used in a multinomial logistic regression to generate a categorical 
#' exposure variable for that item. Each matrix of betas has R rows and number 
#' of columns equal to the number of columns in the design matrix.
#' @seealso [simulate_pop()]
#' @export
#' @examples
#' ## Example 1: X_i ~ C_i
#' # Number of items, exposure levels, latent classes
#' J <- 30; R <- 4; K <- 3
#' # Formula specifying that X_i depends on C_i
#' formula_x <- "~ c_all"
#' # Dataframe with unique values of C_i
#' V_unique <- data.frame(c_all = as.factor(1:K))
#' # Matrix of pattern profiles for each latent class
#' profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J), 
#'                                         rep(3, times = 0.5 * J)),
#'                                  C2 = c(rep(4, times = 0.2 * J), 
#'                                         rep(2, times = 0.8 * J)),
#'                                  C3 = c(rep(3, times = 0.3 * J), 
#'                                         rep(4, times = 0.4 * J),
#'                                         rep(1, times = 0.3 * J))))
#' # True level probability
#' modal_theta_prob <- 0.85
#' # Get matrix of betas for generating C_i
#' beta_list_x <- get_betas_x(profiles = profiles, R = R, 
#'                            modal_theta_prob = modal_theta_prob, 
#'                            formula_x = formula_x, V_unique = V_unique)
#' # Beta matrix for item j=1
#' beta_list_x[[1]]     
#'                       
#' ## Example 2: X_i ~ C_i + S_i
#' # Update formula_x
#' formula_x <- "~ c_all + s_all"
#' # Update beta_list_x by adding in coefficients for S_i to each j matrix
#' beta_list_x <- lapply(1:J, function(j) cbind(beta_list_x[[j]], 
#'                                              s_all = c(0, 0.5, 0, 0)))
#' beta_list_x[[1]]  
#'                      
get_betas_x <- function(profiles, R, modal_theta_prob = 0.85, formula_x, V_unique) {
  # Check errors
  var_terms <- labels(stats::terms(stats::as.formula(formula_x)))
  if (!any(grepl("c_all", var_terms))) {
    stop("formula_x must contain 'c_all' as a covariate so that the exposure is generated from the latent classes")
  } else if (!is.factor(V_unique[["c_all"]])) {
    stop("V_unique must contain the unique values of factor variable c_all")
  }
  if (any(grepl(":", var_terms))) {
    stop("please do not specify interactions in formula_x")
  }
  if (!all(unique(c(profiles)) %in% 1:R)) {
    stop("The item levels in profiles must range from 1 to R, where R is the number of levels for all items")
  }
  
  # Create design matrix with unique values
  design_mat_unique <- stats::model.matrix(stats::as.formula(formula_x), 
                                           data = V_unique)
  
  # Get dimensions and initialize values
  K <- ncol(profiles)
  J <- nrow(profiles)
  Q_full <- ncol(design_mat_unique)
  non_mode <- (1 - modal_theta_prob) / (R - 1)
  mode_div_non <- log(modal_theta_prob / non_mode)
  non_div_mode <- log(non_mode / modal_theta_prob)
  beta_list_x <- vector(mode = "list", length = J)
  
  # For each exposure item j, get beta matrix
  for (j in 1:J) {
    beta_mat_j <- matrix(0, nrow = R, ncol = Q_full)
    colnames(beta_mat_j) <- colnames(design_mat_unique)
    rownames(beta_mat_j) <- 1:R
    theta_j <- profiles[j, ]  # profiles for item j
    # Case when mode for C_i=1 is category 1
    if (theta_j[1] == 1) {  
      # C_i = 1
      beta_mat_j[-1, 1] <- non_div_mode
      # C_i = 2, 3, ...
      for (k in 2:K) {
        # If mode is same as C_i=1 mode, then all 0's
        # Else if mode is different from C_i = 1
        if (theta_j[k] != theta_j[1]) { 
          # beta_rk = M + M*I(r is mode)
          beta_mat_j[-1, k] = mode_div_non
          beta_mat_j[theta_j[k], k] = beta_mat_j[theta_j[k], k] + mode_div_non
        }
      }
    # Case when mode for C_i=1 is not category 1
    } else {  
      # C_i = 1
      beta_mat_j[theta_j[1], 1] <- mode_div_non 
      # C_i = 2, 3, ...
      for (k in 2:K) {
        if (theta_j[k] != theta_j[1]) { # Mode is different from C_i = 1
          if (theta_j[k] == 1) {  # Mode is category 1
            # beta_rk = N + N*I(r is C_i=1 mode)
            beta_mat_j[-1, k] = non_div_mode
            beta_mat_j[theta_j[1], k] = beta_mat_j[theta_j[1], k] + non_div_mode
          } else {  # Mode is not category 1
            # beta_rk = N*I(r is C_i=1 mode) + M*I(r is mode)
            beta_mat_j[theta_j[1], k] = non_div_mode
            beta_mat_j[theta_j[k], k] = mode_div_non
          }
        }
      }
    }
    
    # Add beta matrix to list of betas for all items
    beta_list_x[[j]] <- beta_mat_j
  }
  # Return list of beta matrices
  return(beta_list_x)
  
}


#' Create a categorical variable using multinomial logistic regression
#' 
#' @param beta_mat Coefficient parameters for the linear predictor terms of a 
#' multinomial logistic regression. Kx(Q_full), where K is the number of 
#' categories and Q_full is the number of covariate terms in the regression
#' @param design_mat Regression design matrix. Nx(Q_full)
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
#'   \item{\code{pi_split}}{If `split_dim` is specified, SxK matrix of category 
#'   probabilities split by the levels in `split_dim`, where S is the number of 
#'   levels in the `split_dim` variable, and K is the number of categories as 
#'   specified by `beta_mat`. Otherwise, `NULL` if `split_dim` is `NULL`.}
#' }
#' @seealso [simulate_pop()]
#' @keywords internal
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
#' # (0.3, 0.5, 0.2) for S_i=1 and (0.1, 0.6, 0.3) for S_i=2, giving overall 
#' # probabilities (0.253, 0.522, 0.225)
#' beta_mat_c <- matrix(c(0, 0, 
#'                        0.5, 1.3,
#'                        -0.4, 1.5), byrow = TRUE, nrow = 3, ncol = 2)
#' # Create categorical variable                      
#' out_vars <- create_categ_var(beta_mat = beta_mat_c, design_mat = design_mat_c,
#'                              split_dim = "s_all", V = V)
#' out_vars$categ_var
#' 
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
    lin_preds[, k] <- design_mat %*% as.matrix(beta_mat[k, ]) 
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
#' @param S Number of subpopulations. Default is 2.
#' @param J Number of exposure items. Default is 30.
#' @param K Number of latent classes. Default is 3.
#' @param R Number of exposure categories for all items. Default is 4.
#' @param N_s Population size for each level of categorical variable S_i. 
#' Default is `c(60000, 20000)`, corresponding to a binary S_i with S=2 levels.
#' @param formula_c String specifying formula for multinomial logistic 
#' regression to create category latent class assignment C_i. Default is `"~ s_all"`, 
#' which generates C_i dependent on stratum variable S_i. All variables in the 
#' formula must be `"c_all"`, `"s_all"`, or specified in `V_additional`.
#' @param formula_x String specifying formula for multinomial logistic regression 
#' to create multivariate categorical exposure X_i. Default is `"~ c_all"`, 
#' which generates X_i dependent on latent class C_i. All variables in 
#' the formula must be `"c_all"`, `"s_all"`, or specified in `V_additional`.
#' @param formula_y String specifying formula for logistic regression to create 
#' binary outcome Y_i. Default is `"~ c_all * s_all"`, which generates Y_i dependent 
#' on latent class C_i, stratum S_i, and their interaction. All variables in the 
#' formula must be `"c_all"`, `"s_all"`, or specified in `V_additional`.
#' @param modal_theta_prob Probability between 0 and 1 for the most likely 
#' exposure category, assumed to be the same for all items. For all other
#' categories, 1 - `modal_theta_prob` is evenly split across them. Default is 
#' 0.85, so if there are R=4 exposure categories, then the 3 non-modal 
#' categories occur with probability \eqn{(1-0.85)/3 = 0.05}
#' @param beta_mat_c Coefficient parameters for a multinomial logistic 
#' regression to generate latent class assignment C_i. Kx(Q_full), where K is the 
#' number of categories and Q_full is the number of covariate terms in the regression. 
#' Default is `NULL` and default values are used. See details.
#' @param beta_list_x Coefficient parameters for a multinomial logistic regression 
#' to generate multivariate exposure X_i. List of J matrices, each of dimension 
#' Kx(Q_full). Default is `NULL` and default values are used. See details.
#' @param beta_vec_y Coefficient parameters for a logistic regression to 
#' generate binary outcome Y_i. (Q_full)x1. Default is `NULL` and 
#' default values are used. See details.
#' @param xi_mat_y Alternative to `xi_mat_y` but in mixture reference coding. 
#' Default is `NULL`. If specified, must have K rows and Q columns, where Q is 
#' the number of covariate terms in the regression, excluding all terms 
#' involving C_i.
#' @param V_additional Dataframe with any additional variables other than C_i and
#' S_i to be used to generate variables. Number of rows should be N. Default is `NULL`.
#' @param cluster_size Size of clusters for clustering in the outcome, Y_i, assumed
#' to be the same for all clusters. Default is 80, corresponding to 1000 
#' clusters if N is 80000, assuming exchangeable latent correlation matrix with 
#' correlation 0.5 on the off-diagonals.
#' @param pop_seed Numeric seed for population generation. Default is 1.
#' @param save_res Boolean specifying if results should be saved. Default = `TRUE`.
#' @param save_path String specifying directory and file name to save results,
#' e.g., "~/Documents/stratified_classes". Default is `NULL`.
#' 
#' @details
#' If `NULL` (default) is used for `beta_mat_c`, `beta_list_x`, or `beta_mat_y`, 
#' the following default values are used, corresponding to the scenario with K=3 
#' latent class and S=2 levels for variable S_i: `beta_mat_c` is a 3x2 matrix with
#' values `c(0, 0, 0.5, 1.3, -0.4, 1.5)` by row; `beta_list_x` is a list of J=30
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
#'   \item{\code{S}}{Number of stratum (i.e., levels of S_i)}
#'   \item{\code{N_s}}{Vector of population sizes for levels of S_i. Sx1}
#'   \item{\code{true_K}}{True number of latent classes, K}
#'   \item{\code{true_Ai}}{`NULL` or vector of additional continuous variable if 
#'   `a_all` is provided in `V_additional`. Nx1}
#'   \item{\code{true_Bi}}{`NULL` or vector of additional binary variable if 
#'   `b_all` is provided in `V_additional`. Nx1}
#'   \item{\code{true_Si}}{Vector of true stratum indicators. Nx1}
#'   \item{\code{true_Ci}}{Vector of true latent class indicators. Nx1}
#'   \item{\code{true_pi}}{Vector of true pi values overall in the population. Kx1}
#'   \item{\code{true_pi_s}}{`NULL` or SxK matrix of true pi values within each 
#'   level of S_i}
#'   \item{\code{X_data}}{Matrix of multivariate categorical exposure for all 
#'   individuals. NxJ}
#'   \item{\code{true_global_patterns}}{Matrix of true global exposure patterns 
#'   defined by modal category. JxK}
#'   \item{\code{true_global_thetas}}{Array of true thetas. JxKxR}
#'   \item{\code{Y_data}}{Vector of binary outcome for all individuals. Nx1}
#'   \item{\code{cluster_id}}{Vector of true cluster indicators. Nx1}
#'   \item{\code{cluster_size}}{Cluster size}
#'   \item{\code{true_xi}}{Matrix of probit regression coefficients in mixture 
#'   reference coding. KxQ, where Q is the number of covariate terms in the 
#'   regression, excluding all terms involving C_i}
#'   \item{\code{true_Phi}}{Vector of true outcome probabilities for all 
#'   individuals. Nx1}
#'   \item{\code{true_Phi_mat}}{`NULL` or matrix of true outcome probabilities
#'   for individuals aggregated by C_i and S_i. KxS}
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
#' \dontrun{
#' ### Example 1: Default values
#' sim_pop <- simulate_pop(save_res = FALSE)
#' 
#' ### Example 2: Similar to default but smaller population size
#' # Population size and strata dimensions
#' N = 800; S = 2; N_s = c(600, 200)
#' 
#' # Generate C_i ~ S_i
#' K <- 3  
#' formula_c <- "~ s_all"
#' V_unique <- data.frame(s_all = as.factor(1:S))
#' pi_mat <- matrix(c(0.3, 0.5, 0.2,   
#'                    0.1, 0.6, 0.3), 
#'                  byrow = TRUE, nrow = S, ncol = K)
#' beta_mat_c <- get_betas_c(pi_mat = pi_mat, formula_c = formula_c, 
#'                           V_unique = V_unique)
#' 
#' # Generate X_i ~ C_i
#' J <- 30; R <- 4
#' formula_x <- "~ c_all"
#' V_unique <- data.frame(c_all = as.factor(1:K))
#' profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J), 
#'                                         rep(3, times = 0.5 * J)),
#'                                  C2 = c(rep(4, times = 0.2 * J), 
#'                                         rep(2, times = 0.8 * J)),
#'                                  C3 = c(rep(3, times = 0.3 * J), 
#'                                         rep(4, times = 0.4 * J),
#'                                         rep(1, times = 0.3 * J))))
#' modal_theta_prob <- 0.85
#' beta_list_x <- get_betas_x(profiles = profiles, R = R, 
#'                            modal_theta_prob = modal_theta_prob, 
#'                            formula_x = formula_x, V_unique = V_unique)
#' 
#' # Generate Y_i ~ C_i + S_i + C_i:S_i
#' formula_y <- "~ c_all * s_all"
#' beta_vec_y <- c(1, -0.7, -1.5, -0.5, -0.5, -0.3)
#' cluster_size <- 80
#' 
#' # Simulate population
#' pop_seed <- 1  # Set seed
#' sim_pop <- simulate_pop(N = N, J = J, K = K, R = R, N_s = N_s,
#'                          modal_theta_prob = modal_theta_prob, 
#'                          formula_c = formula_c, formula_x = formula_x, 
#'                          formula_y = formula_y, beta_mat_c = beta_mat_c, 
#'                          beta_list_x = beta_list_x, beta_vec_y = beta_vec_y, 
#'                          cluster_size = cluster_size, 
#'                          pop_seed = pop_seed, save_res = FALSE)    
#'                                               
#' ### Example 2: Selection-dependent pattern profiles     
#' # Generate X_i ~ C_i + S_i
#'  formula_x <- "~ c_all + s_all"
#'  beta_list_x <- lapply(1:J, function(j) cbind(beta_list_x[[j]], 
#'                                               s_all = c(0, 0.5, 0, 0)))
#' 
#' # Create population
#' sim_pop <- simulate_pop(N = N, S=S, J = J, K = K, R = R, N_s = N_s,
#'                         modal_theta_prob = modal_theta_prob, 
#'                         formula_c = formula_c, formula_x = formula_x, 
#'                         formula_y = formula_y, beta_mat_c = beta_mat_c, 
#'                         beta_list_x = beta_list_x, beta_vec_y = beta_vec_y, 
#'                         cluster_size = cluster_size, 
#'                         pop_seed = pop_seed, save_res = FALSE)
#' 
#' ### Example 3: Additional effect modifiers for Y_i and no clustering      
#' # Continuous variable A_i for age centered about 0
#' a_all <- stats::rnorm(n = N, mean = 0, sd = 5)
#' # Binary variable B_i for physically inactive or active
#' b_all <- as.factor(stats::rbinom(n = N, size = 1, prob = 0.3) + 1)
#' # Create dataframe of additional variables A_i and B_i
#' V_additional <- data.frame(a_all, b_all)    
#'               
#' # Generate Y_i ~ C_i + S_i + A_i + B_i + C_i:S_i + C_i:A_i + C_i:B_i
#' formula_y <- "~ c_all * (s_all + a_all + b_all)"
#' beta_vec_y <- c(1, -0.7, -1.5, -0.5, -0.5, -0.3, -0.04, 0.09, 0.08, 0.4, 
#'                 -0.7, -0.6)
#' cluster_size <- 1
#'                 
#' # Create population
#' sim_pop <- simulate_pop(N = N, S=S, J = J, K = K, R = R, N_s = N_s,
#'                         modal_theta_prob = modal_theta_prob, 
#'                         formula_c = formula_c, formula_x = formula_x, 
#'                         formula_y = formula_y, beta_mat_c = beta_mat_c, 
#'                         beta_list_x = beta_list_x, beta_vec_y = beta_vec_y, 
#'                         V_additional = V_additional, cluster_size = cluster_size, 
#'                         pop_seed = pop_seed, save_res = FALSE)
#' }
simulate_pop <- function(N = 80000, S = 2, J = 30, K = 3, R = 4, 
                         N_s = c(60000, 20000),  modal_theta_prob = 0.85, 
                         formula_c = "~ s_all", 
                         formula_x = "~ c_all", 
                         formula_y = "~ c_all * s_all", 
                         beta_mat_c = NULL, beta_list_x = NULL,
                         beta_vec_y = NULL, xi_mat_y = NULL, 
                         V_additional = NULL, 
                         cluster_size = 80, pop_seed = 1, 
                         save_res = TRUE, save_path = NULL) {

  # Set seed
  set.seed(pop_seed)
  #================== Check errors =============================================
  
  # Catch errors
  # Check number of strata
  if (length(N_s) != S) {
    stop("N_s must be a numeric vector of length S")
  }
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
  non_mode <- (1 - clust_mode) / (R - 1)
  
  #================ Create S_i variable ==========================================
  # Create strata
  s_all <- unlist(sapply(1:S, function(x) rep(x, times = N_s[x])))
  # Dataframe of strata and additional variables
  if (!is.null(V_additional)) {
    V <- data.frame(s_all = as.factor(s_all), V_additional)
  } else {
    V <- data.frame(s_all = as.factor(s_all))
  }
  
  #================ Create latent class assignments C_i variable =================
  # Set defaults for generating latent class dependent on S_i and check beta_mat_c
  # for errors
  if (is.null(beta_mat_c)) {
    if (K != 3 | S != 2) {
      stop("If default for beta_mat_c is to be used, K must be equal to 3 and S 
           must be equal to 2.")
    }
    # Corresponds to an overall true_pi ~= (0.253, 0.522, 0.225)
    # Matrix of class membership probabilities for each level of S_i
    pi_mat <- matrix(c(0.3, 0.5, 0.2,   # class membership probs for S_i=1
                       0.1, 0.6, 0.3),  # class membership probs for S_i=2
                     byrow = TRUE, nrow = S, ncol = K)
    beta_mat_c <- get_betas_c(pi_mat = pi_mat, formula_c = formula_c, 
                              V_unique = V)
  } else if (!is.matrix(beta_mat_c)) {
    # Check that beta_mat_c is a matrix
    stop("beta_mat_c must be a matrix")
  }
  # Check that all variables in formula_c are available, removing interactions
  regr_vars <- labels(stats::terms(stats::as.formula(formula_c)))
  regr_vars <- regr_vars[!(stringr::str_detect(regr_vars, ":"))] 
  if (!all(regr_vars %in% colnames(V))) {
    stop("all variables in formula_c other than s_all must be provided in V_additional")
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
  
  #================ Create multivariate exposure X_i variable ====================
  # Add latent class to data frame of variables
  V <- data.frame(c_all = as.factor(c_all), V)
  # Get design matrix
  design_mat_x <- stats::model.matrix(stats::as.formula(formula_x), data = V)
  
  # Set defaults for generating exposure X_i dependent on C_i and S_i and check
  # beta_list_x for errors
  if (is.null(beta_list_x)) {
    if (K != 3) {
      stop("If default for beta_list_x is to be used, K must be equal to 3.")
    }
    profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                          rep(3, times = 0.5 * J)),
                                   C2 = c(rep(4, times = 0.2 * J),
                                          rep(2, times = 0.8 * J)),
                                   C3 = c(rep(3, times = 0.3 * J),
                                          rep(4, times = 0.4 * J),
                                          rep(1, times = 0.3 * J))))
    beta_list_x <- get_betas_x(profiles = profiles, R = R, 
                              modal_theta_prob = modal_theta_prob, 
                              formula_x = formula_x, V_unique = V)
  } else if (!is.list(beta_list_x) | !is.matrix(beta_list_x[[1]])) {
    # Check that beta_list_x is a list of matrices
    stop("beta_list_x must be a list of matrices")
  } else if (ncol(beta_list_x[[1]]) != ncol(design_mat_x)) {
    stop(paste0("Each matrix in beta_list_x must have columns corresponding to '", 
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
    out_vars_x <- create_categ_var(beta_mat = beta_list_x[[j]], 
                                   design_mat = design_mat_x, 
                                   split_dim = "c_all", V = V)
    X_data[, j] <- out_vars_x$categ_var
    true_global_thetas[j, , ] <- out_vars_x$pi_split
    true_global_patterns[j, ] <- apply(out_vars_x$pi_split, 1, which.max)
  }
  
  #================ Create binary Y_i outcome variable ===========================
  # Set defaults for generating exposure Y_i dependent on C_i and S_i and check
  # xi_mat_y and beta_vec_y for errors
  if (is.null(xi_mat_y)) {
    if (is.null(beta_vec_y)) {
      if (K != 3 | S != 2) {
        stop("If default for beta_vec_y is to be used, K must be equal to 3 and S 
           must be equal to 2.")
      }
      # Corresponds to 3x2 xi_mat_y with c(1, -0.5, 
      #                                    0.3, -1, 
      #                                    -0.5, -0.8) by row
      # and probabilities 0.84, 0.62, 0.31 for S_i=1 and 0.69, 0.24, 0.1 for S_i=2
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
    # xi1*I(C_i=1) + xi2*I(C_i=1,S_i=2) +
    # + xi3*I(C_i=2) + xi4*I(C_i=2,S_i=2)
    # + xi5*I(C_i=3) + xi6*I(C_i=3,S_i=2)
    # Number of covariate terms in the regression, excluding all terms with C_i
    cov_terms <- colnames(stats::model.matrix(stats::as.formula(formula_y)))
    cov_terms_no_c <- cov_terms[!stringr::str_detect(cov_terms, "c_all")]
    Q <- length(cov_terms_no_c)
    xi_mat_y <- convert_ref_to_mix(K = K, Q = Q, est_beta = beta_vec_y)$est_xi
  }
  
  # Get vector of individual linear predictors
  lin_pred <- design_mat_y %*% beta_vec_y
  # Get vector of individual underlying outcome probabilities
  true_Phi_under <- stats::pnorm(lin_pred)
  
  ## Generate outcome Y_i depending on clustering
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
  # output outcome probabilities individually and aggregated by C_i and S_i
  y_covs <- labels(stats::terms(stats::as.formula(formula_y)))
  if (all(y_covs %in% c("c_all", "s_all", "c_all:s_all"))) {
    pop_inds <- 1:N
    # Matrix of finite population outcome probs aggregated by C_i and S_i
    true_Phi_mat <- matrix(NA, nrow=K, ncol=S)
    # Observed finite population probs for each individual based on C_i and S_i
    true_Phi <- numeric(N)
    for (h in 1:S) {
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
  
  sim_pop <- list(N = N, J = J, R = R, S = S, N_s = N_s, true_K = K,
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
#' @param strat Boolean specifying whether to perform stratification by S_i.
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
#'   \item{\code{S}}{Number of stratum (i.e., levels of S_i)}
#'   \item{\code{N_s}}{Vector of population sizes for levels of S_i. Sx1}
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
#'   \item{\code{true_pi_s}}{`NULL` or SxK matrix of true pi values within each 
#'   level of S_i}
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
#'   reference coding. KxQ, where Q is the number of covariate terms in the 
#'   regression, excluding all terms involving C_i}
#'   \item{\code{true_Phi}}{Vector of true outcome probabilities for sampled 
#'   individuals. nx1}
#'   \item{\code{true_Phi_mat}}{`NULL` or matrix of true outcome probabilities
#'   for individuals aggregated by C_i and S_i. KxS}
#' }
#' If `save_res = TRUE` (default), also saves `sim_samp` as 
#' `[save_path]_sim_samp.RData`. 
#' 
#' @seealso [simulate_pop()] 
#' @export
#' @examples 
#' \dontrun{
#' ### Example 1: Default values
#' # Create population
#' sim_pop <- simulate_pop(save_res = FALSE)
#' sim_samp <- simulate_samp(sim_pop = sim_pop, samp_prop = 0.05, strat = TRUE,
#'                           strat_dist = c(0.5, 0.5), clust = TRUE,
#'                           samp_seed = 101, save_res = FALSE)
#' }
simulate_samp <- function(sim_pop, samp_prop = 0.05, samp_size = NULL, 
                          strat = TRUE, strat_dist = c(0.5, 0.5), clust = TRUE, 
                          samp_seed = 101, save_res = TRUE, save_path) {
  # Set seed
  set.seed(samp_seed)
  #================ Check errors ===============================================
  # Check clustering
  if (sim_pop$cluster_size == 1 & clust) {
    warning("cluster sampling is requested but population does not contain 
                clustered data")
  }
  # Check stratification distribution
  if (strat) {
    if (is.null(strat_dist)) {
      stop("strat_dist must be specified if strat is TRUE")
    } else if (sum(strat_dist) != 1) {
      stop("strat_dist must be a vector of length S whose elements sum to 1")
    }
  }
  # Check save_res and save_path
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
             file named 'stratified_classes_sim_samp.RData'")
      }
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
    if (length(strat_dist) != sim_pop$S) {
      stop(paste0("strat_dist must be a vector of length S = ", sim_pop$S, 
                  ", to match the number of levels in stratifying variable S_i"))
    }
    # Get stratum sample sizes
    n_s <- samp_size * strat_dist
    for (h in 1:sim_pop$S) {
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
      if (n_clus %% 1 != 0) {
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
    for (h in 1:sim_pop$S) {
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
                   N = sim_pop$N, J = sim_pop$J, R = sim_pop$R, S = sim_pop$S, 
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


