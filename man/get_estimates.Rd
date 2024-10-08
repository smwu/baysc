% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_functions.R
\name{get_estimates}
\alias{get_estimates}
\title{Get posterior estimates}
\usage{
get_estimates(MCMC_out, post_MCMC_out, n, J, V, y_all, x_mat)
}
\arguments{
\item{MCMC_out}{Output from \code{\link[=run_MCMC_Rcpp]{run_MCMC_Rcpp()}} containing:
\describe{
\item{\code{pi_MCMC}}{Matrix of posterior samples for pi. (n_iter)xK}
\item{\code{theta_MCMC}}{Array of posterior samples for theta. (n_iter)xJxKxR}
\item{\code{xi_MCMC}}{Array of posterior samples for xi. (n_iter)xKxQ}
\item{\code{c_all_MCMC}}{Matrix of posterior samples for c_all. (n_iter)xn}
\item{\code{z_all_MCMC}}{Matrix of posterior samples for z_all. (n_iter)xn}
\item{\code{loglik_MCMC}}{Vector of posterior samples for log-likelihood. (n_iter)x1}
}}

\item{post_MCMC_out}{output from \code{\link[=post_process]{post_process()}} containing:
\describe{
\item{\code{K_med}}{Median, across iterations, of number of classes with at least \code{class_cutoff} percent of individuals}
\item{\code{pi}}{Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)}
\item{\code{theta}}{Array of reduced and relabeled posterior samples for theta. (n_iter)xJx(K_med)xR}
\item{\code{xi}}{Array of reduced and relabeled posterior samples for xi. (n_iter)x(K_med)xQ}
\item{\code{dendrogram}}{Hierarchical clustering dendrogram used for relabeling}
}}

\item{n}{Number of individuals}

\item{J}{Number of exposure items}

\item{V}{Regression design matrix without class assignment. nxQ}

\item{y_all}{Vector of outcomes. nx1}

\item{x_mat}{Matrix of multivariate categorical exposures. nxJ}
}
\value{
Returns list \code{estimates} containing:
\describe{
\item{\code{K_red}}{Number of unique classes}
\item{\code{pi_red}}{Matrix of final posterior samples for pi. Mx(K_red),
where M is the number of MCMC iterations after burn-in and thinning.}
\item{\code{theta_red}}{Array of final posterior samples for theta. MxJx(K_red)xR}
\item{\code{xi_red}}{Array of final posterior samples for xi. Mx(K_red)xQ}
\item{\code{pi_med}}{Vector of posterior median estimates for pi. (K_red)x1}
\item{\code{theta_med}}{Array of posterior median estimates for theta. Jx(K_red)xR}
\item{\code{xi_med}}{Matrix of posterior median estimates for xi. (K_red)xQ}
\item{\code{Phi_med}}{Vector of final individual outcome probabilities. nx1}
\item{\code{c_all}}{Vector of final individual class assignments. nx1}
\item{\code{pred_class_probs}}{Matrix of individual posterior class probabilities. nx(K_red)}
\item{\code{loglik_med}}{Vector of final indiviudal log-likehoods. nx1}
}
}
\description{
\code{get_estimates} obtains posterior parameter samples and estimates prior to
variance adjustment
}
\details{
First, duplicate classes that have the same modal exposure categories
for all items are combined to obtain the number of unique classes, \code{K_red}.
Parameters are then renormalized for the unique classes and posterior median
estimates are computed across MCMC iterations. Using these median estimates,
class assignments \code{c_all}, the regression mean, and the individual
log-likelihood are derived.
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

# Then run post-process relabeling
post_MCMC_out <- post_process(MCMC_out = MCMC_out, J = J, R = R, Q = Q,
class_cutoff = 0.05)

# Then obtain posterior estimates
estimates <- get_estimates(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out,
                           n = n, J = J, V = V, y_all = y_all, x_mat = x_mat)

}
\references{
Williams, M. R. and Savitsky, T. D. (2021). Uncertainty estimation for
pseudo-bayesian inference under complex sampling. International Statistical
Review 89, 72–107.
}
\seealso{
\code{\link[=run_MCMC_Rcpp]{run_MCMC_Rcpp()}} \code{\link[=post_process]{post_process()}} \code{\link[=swolca_var_adjust]{swolca_var_adjust()}} \code{\link[=swolca]{swolca()}}
}
