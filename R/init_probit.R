#' Initialize probit model
#'
#' @description
#' `init_probit` initializes priors and variables for the probit regression model 
#' given hyperparameters.
#'
#' @inheritParams init_OLCA
#' @param q Number of regression covariates excluding class assignment
#' @param V Regression design matrix without class assignment. nxq
#' @param mu0 List of K qx1 vectors of mean hyperparameters for regression 
#' coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}. 
#' @param Sig0 List of K qxq matrices of variance hyperparameters for regression 
#' coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}. 
#' @param y_all Vector of outcomes. nx1
#' @param c_all Vector of random initial class assignments. nx1
#' 
#' @details
#' First, regression coefficients \eqn{\xi} are initialized by drawing independently
#' for each latent class from a Multivariate Normal distribution with mean vector 
#' hyperparameter \eqn{\mu_0} drawn from a Normal(0,1) hyperprior and variance 
#' matrix hyperparameter \eqn{\Sigma_0} set to be a diagonal matrix with 
#' diagonal components drawn from InvGamma(shape=5/2, scale=2/5) distributions.
#' Then, the latent probit variable `z_all` is initialized by drawing from a 
#' Truncated Normal distribution with mean \eqn{V\xi} and variance 1 and 
#' truncation boundary 0, where negative values are drawn if the outcome is 0 and
#' positive values are drawn if the outcome is 1. 
#'
#' @return
#' Returns list `probit_params` containing:
#' \describe{
#'   \item{\code{xi}}{Matrix parameter xi for probit regression coefficients. Kxq}
#'   \item{\code{z_all}}{Vector of latent variables in the probit model. nx1}
#' }
#'
#' @seealso [init_OLCA()] [swolca()] [solca()]
#' 
#' @importFrom LaplacesDemon rmvn rinvgamma
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats qnorm rnorm
#' @export
#'
#' @examples
#' # Load data and obtain relevant variables
#' data("sim_data")
#' data_vars <- sim_data
#' x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
#' y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
#' n <- dim(x_mat)[1]        # Number of individuals
#' J <- dim(x_mat)[2]        # Number of exposure items
#' R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
#'              function(x) length(unique(x)))  
#' R <- max(R_j)             # Maximum number of exposure categories across items
#' 
#' # Probit model only includes latent class
#' V_data <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
#' glm_form <- "~ 1"
#' # Obtain probit regression design matrix without class assignment
#' V <- model.matrix(as.formula(glm_form), data = V_data)
#' # Number of regression covariates excluding class assignment
#' q <- ncol(V)  
#' 
#' # Set hyperparameters
#' K <- 30
#' alpha <- rep(1, K) / K
#' eta <- matrix(0.01, nrow = J, ncol = R) 
#' for (j in 1:J) {
#'   eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
#' }
#' mu0 <- Sig0 <- vector("list", K)
#' for (k in 1:K) {
#'   # MVN(0,1) hyperprior for prior mean of xi
#'   mu0[[k]] <- stats::rnorm(n = q)
#'   # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
#'   # components and mean variance 2.5 for a weakly informative prior on xi
#'   Sig0[[k]] <- diag(LaplacesDemon::rinvgamma(n = q, shape = 3.5, scale = 6.25), 
#'   nrow = q, ncol = q)
#' }
#' 
#' # First initialize OLCA params
#' OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)
#' 
#' # Initialize probit model
#' probit_params <- init_probit(K = K, n = n, q = q, V = V, mu0 = mu0, 
#' Sig0 = Sig0, y_all = y_all, c_all = OLCA_params$c_all)
#' # probit_params
#' 
init_probit <- function(K, n, q, V, mu0, Sig0, y_all, c_all) {
  # Initialize variables
  xi <- matrix(NA, nrow = K, ncol = q)
  z_all <- lin_pred <- numeric(n)
  
  # Prior for xi. Same prior for each class
  for (k in 1:K) {
    xi[k, ] <- LaplacesDemon::rmvn(n = 1, mu = mu0[[k]], Sigma = Sig0[[k]])
  }
  
  # Initialize probit model latent variable, z, following Albert and Chib (1993)
  # Linear predictor using covariate values and class assignment for each individual
  for (i in 1:n) {
    lin_pred[i] <- V[i, ] %*% xi[c_all[i], ]  # V[i]*xi[c_i]
  }
  # For cases, z_i ~ TruncNormal(low=0, high=Inf, mean=V*xi[c_i], var=1)
  z_all[y_all == 1] <- truncnorm::rtruncnorm(n = 1, a = 0, b = Inf, 
                                  mean = lin_pred[y_all == 1], sd = 1)
  # For controls, z_i ~ TruncNormal(low=-Inf, high=0, mean=V*xi[c_i], var=1)
  z_all[y_all == 0] <- truncnorm::rtruncnorm(n = 1, a = -Inf, b = 0, 
                                  mean = lin_pred[y_all == 0], sd = 1)
  # Control extreme values
  z_all[z_all == Inf] <- stats::qnorm(1 - (1e-10))
  z_all[z_all == -Inf] <- stats::qnorm(1e-10)
  
  probit_params <- list(xi = xi, z_all = z_all)
  return(probit_params)
}