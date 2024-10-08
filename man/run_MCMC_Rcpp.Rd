% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_functions.R
\name{run_MCMC_Rcpp}
\alias{run_MCMC_Rcpp}
\title{Run MCMC to get posterior samples}
\usage{
run_MCMC_Rcpp(
  OLCA_params,
  probit_params,
  n_runs,
  burn,
  thin,
  K,
  J,
  R,
  n,
  Q,
  w_all,
  x_mat,
  y_all,
  V,
  alpha,
  eta,
  mu0,
  Sig0,
  update = 10
)
}
\arguments{
\item{OLCA_params}{Output list from \code{\link[=init_OLCA]{init_OLCA()}} containing:
\describe{
\item{\code{pi}}{Vector parameter pi for class membership probabilities. Kx1}
\item{\code{theta}}{Array parameter theta for item level probabilities. JxKxR}
\item{\code{c_all}}{Vector of random initial class assignments. nx1}
}}

\item{probit_params}{Output list from \code{\link[=init_probit]{init_probit()}} containing:
\describe{
\item{\code{xi}}{Matrix parameter xi for probit regression coefficients. KxQ}
\item{\code{z_all}}{Vector of latent variables in the probit model. nx1}
}}

\item{n_runs}{Number of MCMC iterations. Default is 20000.}

\item{burn}{Number of MCMC iterations to drop as a burn-in period. Default is 10000.}

\item{thin}{Thinning factor for MCMC iterations. Default is 5.}

\item{K}{Number of classes}

\item{J}{Number of exposure items}

\item{R}{Maximum number of exposure categories}

\item{n}{Number of individuals}

\item{Q}{Number of regression covariates excluding class assignment}

\item{w_all}{Weights normalized to sum to n. nx1}

\item{x_mat}{Matrix of multivariate categorical exposures. nxJ}

\item{y_all}{Vector of outcomes. nx1}

\item{V}{Regression design matrix without class assignment. nxQ}

\item{alpha}{Kx1 vector of hyperparameters for prior for class membership
probabilities \eqn{\pi}}

\item{eta}{JxR matrix of hyperparameters for prior for item consumption
level probabilities \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class
\eqn{k}, assumed to be the same across classes.}

\item{mu0}{List of K Qx1 vectors of mean hyperparameters for regression
coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}.}

\item{Sig0}{List of K QxQ matrices of variance hyperparameters for regression
coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}.}

\item{update}{Number specifying that MCMC progress updates should be printed
every \code{update} iterations. Default is 10.}
}
\value{
Returns list \code{MCMC_out} containing:
\describe{
\item{\code{pi_MCMC}}{Matrix of posterior samples for pi. (n_iter)xK}
\item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xJxKxR}
\item{\code{xi_MCMC}}{Array of posterior samples for xi. (n_iter)xKxQ}
\item{\code{c_all_MCMC}}{Matrix of posterior samples for c_all. (n_iter)xn}
\item{\code{z_all_MCMC}}{Matrix of posterior samples for z_all. (n_iter)xn}
\item{\code{loglik_MCMC}}{Vector of posterior samples for log-likelihood. (n_iter)x1}
}
}
\description{
\code{run_MCMC_Rcpp} runs the Gibbs sampler MCMC algorithm using Rcpp to obtain
posterior samples.
}
\details{
A Gibbs sampler updates the parameters and variables in the following order:
\eqn{\pi}, \code{c_all}, \eqn{\theta}, \eqn{\xi}, \code{z_all}. Class assignments
are permuted every 10 iterations to encourage mixing, according to a random
permutation sampler (Fruhwirth-Schnatter, 2001).
}
\examples{
   
# Load data and obtain relevant variables
data("sim_data")
data_vars <- sim_data
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
cluster_id <- data_vars$cluster_id  # Cluster indicators, nx1
sampling_wt <- data_vars$sample_wt

# Obtain dimensions
n <- dim(x_mat)[1]        # Number of individuals
J <- dim(x_mat)[2]        # Number of exposure items
R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
             function(x) length(unique(x)))  
R <- max(R_j)             # Maximum number of exposure categories across items
# Obtain normalized weights
kappa <- sum(sampling_wt) / n   # Weights norm. constant
w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1

# Probit model only includes latent class
V_data <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
glm_form <- "~ 1"
# Obtain probit regression design matrix without class assignment
V <- model.matrix(as.formula(glm_form), data = V_data)
# Number of regression covariates excluding class assignment
Q <- ncol(V)  

# Set hyperparameters
K <- 30
alpha <- rep(1, K) / K
eta <- matrix(0.01, nrow = J, ncol = R) 
for (j in 1:J) {
  eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
}
mu0 <- Sig0 <- vector("list", K)
for (k in 1:K) {
  # MVN(0,1) hyperprior for prior mean of xi
  mu0[[k]] <- stats::rnorm(n = Q)
  # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated
  # components and mean variance 2.5 for a weakly informative prior on xi
  Sig0[[k]] <- diag(1/(stats::rgamma(n = Q, shape = 3.5, rate = 6.25)), 
  nrow = Q, ncol = Q)
}

# First initialize OLCA params
OLCA_params <- init_OLCA(K = K, n = n, J = J, R = R, alpha = alpha, eta = eta)

# Then initialize probit params 
probit_params <- init_probit(K = K, n = n, Q = Q, V = V, mu0 = mu0, 
Sig0 = Sig0, y_all = y_all, c_all = OLCA_params$c_all)

# Then run MCMC sampling
MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, 
probit_params = probit_params, n_runs = 50, burn = 25, thin = 5,
K = K, J = J, R = R, n = n, Q = Q, w_all = w_all, x_mat = x_mat, 
y_all = y_all, V = V, alpha = alpha, eta = eta, Sig0 = Sig0, mu0 = mu0)
# MCMC_out

}
\references{
Fruhwirth-Schnatter, S. (2001). Markov chain monte carlo estimation of
classical and dynamic switching and mixture models. Journal of the American
Statistical Association 96, 194–209.
}
\seealso{
\code{\link[=post_process]{post_process()}} \code{\link[=get_estimates]{get_estimates()}} \code{\link[=swolca_var_adjust]{swolca_var_adjust()}} \code{\link[=swolca]{swolca()}}
\code{\link[=run_MCMC_Rcpp_wolca]{run_MCMC_Rcpp_wolca()}}
}
