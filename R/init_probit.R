#' Initialize probit model
#'
#' @description
#' `init_probit` initializes priors and variables for the probit regression model 
#' given hyperparameters.
#'
#' @inheritParams init_OLCA
#' @param q Number of regression covariates excluding class assignment
#' @param V Regression design matrix without class assignment. nxq
#' @param mu0 List of vectors of mean hyperparameters for xi. List of K qx1 vectors
#' @param Sig0 List of matrices of variance hyperparameters for xi. List of K qxq matrices
#' @param y_all Vector of outcomes. nx1
#' @param c_all Vector of random initial class assignments. nx1
#'
#' @return
#' Returns list `probit_params` containing:
#' \describe{
#'   \item{\code{xi}}{Matrix parameter xi for probit regression coefficients. Kxq}
#'   \item{\code{z_all}}{Vector of latent variables in the probit model. nx1}
#' }
#'
#' @export
#'
#' @examples
#' # Probit model only includes latent class 
#' K <- 30; n <- 400; q <- 1
#' V <- matrix(1, nrow = n)
#' mu0 <- Sig0 <- vector("list", K_max)
#' for (k in 1:K_max) {
#'   # MVN(0,1) hyperprior for prior mean of xi
#'   mu0[[k]] <- rnorm(n = q)
#'   # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
#'   # components and mean variance 2.5 for a weakly informative prior on xi
#'   Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
#' }
#' y_all <- rbinom(n = n, size = 1, prob = 0.5)
#' c_all <- c(rep(1, 100), rep(2, 100), rep(3, 100), rep(4, 100))
#' init_probit(K = K, n = n, q = q, V = V, mu0 = mu0, Sig0 = Sig0, 
#' y_all = y_all, c_all = c_all)
#' 
init_probit <- function(K, n, q, V, mu0, Sig0, y_all, c_all) {
  # Initialize variables
  xi <- matrix(NA, nrow = K, ncol = q)
  z_all <- lin_pred <- numeric(n)
  
  # Prior for xi. Same prior for each class
  for (k in 1:K) {
    xi[k, ] <- rmvn(n = 1, mu = mu0[[k]], Sigma = Sig0[[k]])
  }
  
  # Initialize probit model latent variable, z, following Albert and Chib (1993)
  # Linear predictor using covariate values and class assignment for each individual
  for (i in 1:n) {
    lin_pred[i] <- V[i, ] %*% xi[c_all[i], ]  # V[i]*xi[c_i]
  }
  # For cases, z_i ~ TruncNormal(low=0, high=Inf, mean=V*xi[c_i], var=1)
  z_all[y_all == 1] <- rtruncnorm(n = 1, a = 0, b = Inf, 
                                  mean = lin_pred[y_all == 1], sd = 1)
  # For controls, z_i ~ TruncNormal(low=-Inf, high=0, mean=V*xi[c_i], var=1)
  z_all[y_all == 0] <- rtruncnorm(n = 1, a = -Inf, b = 0, 
                                  mean = lin_pred[y_all == 0], sd = 1)
  # Control extreme values
  z_all[z_all == Inf] <- qnorm(1 - (1e-10))
  z_all[z_all == -Inf] <- qnorm(1e-10)
  
  probit_params <- list(xi = xi, z_all = z_all)
  return(probit_params)
}