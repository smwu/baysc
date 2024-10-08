% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_functions.R
\name{init_probit}
\alias{init_probit}
\title{Initialize probit model}
\usage{
init_probit(K, n, Q, V, mu0, Sig0, y_all, c_all)
}
\arguments{
\item{K}{Number of classes}

\item{n}{Number of individuals}

\item{Q}{Number of regression covariates excluding class assignment}

\item{V}{Regression design matrix without class assignment. nxQ}

\item{mu0}{List of K Qx1 vectors of mean hyperparameters for regression
coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}.}

\item{Sig0}{List of K QxQ matrices of variance hyperparameters for regression
coefficients \eqn{\xi_{k\cdot}} for each class \eqn{k}.}

\item{y_all}{Vector of outcomes. nx1}

\item{c_all}{Vector of random initial class assignments. nx1}
}
\value{
Returns list \code{probit_params} containing:
\describe{
\item{\code{xi}}{Matrix parameter xi for probit regression coefficients. KxQ}
\item{\code{z_all}}{Vector of latent variables in the probit model. nx1}
}
}
\description{
\code{init_probit} initializes priors and variables for the probit regression model
given hyperparameters.
}
\details{
First, regression coefficients \eqn{\xi} are initialized by drawing independently
for each latent class from a Multivariate Normal distribution with mean vector
hyperparameter \eqn{\mu_0} drawn from a Normal(0,1) hyperprior and variance
matrix hyperparameter \eqn{\Sigma_0} set to be a diagonal matrix with
diagonal components drawn from InvGamma(shape=5/2, scale=2/5) distributions.
Then, the latent probit variable \code{z_all} is initialized by drawing from a
Truncated Normal distribution with mean \eqn{V\xi} and variance 1 and
truncation boundary 0, where negative values are drawn if the outcome is 0 and
positive values are drawn if the outcome is 1.
}
\examples{
# Load data and obtain relevant variables
data("sim_data")
data_vars <- sim_data
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
n <- dim(x_mat)[1]        # Number of individuals
J <- dim(x_mat)[2]        # Number of exposure items
R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
             function(x) length(unique(x)))  
R <- max(R_j)             # Maximum number of exposure categories across items

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

# Initialize probit model
probit_params <- init_probit(K = K, n = n, Q = Q, V = V, mu0 = mu0, 
Sig0 = Sig0, y_all = y_all, c_all = OLCA_params$c_all)
# probit_params

}
\seealso{
\code{\link[=init_OLCA]{init_OLCA()}} \code{\link[=swolca]{swolca()}}
}
