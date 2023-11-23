% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wolca.R
\name{wolca}
\alias{wolca}
\title{Run the WOLCA model}
\usage{
wolca(
  x_mat,
  y_all,
  sampling_wt,
  cluster_id,
  stratum_id,
  V_data = NULL,
  run_sampler = "both",
  glm_form,
  K_max = 30,
  adapt_seed = NULL,
  class_cutoff = 0.05,
  alpha_adapt = NULL,
  eta_adapt = NULL,
  alpha_fixed = NULL,
  eta_fixed = NULL,
  K_fixed = NULL,
  fixed_seed = NULL,
  n_runs = 20000,
  burn = 10000,
  thin = 5,
  ci_level = 0.95,
  save_res = TRUE,
  save_path = NULL
)
}
\arguments{
\item{x_mat}{Categorical exposure matrix. nxJ}

\item{y_all}{Vector of outcomes. nx1}

\item{sampling_wt}{Vector of survey sampling weights. nx1. If no sampling
weights are available, set this to a vector of ones.}

\item{cluster_id}{Vector of individual cluster IDs. nx1. Default is \code{NULL},
indicating each individual is their own cluster.}

\item{stratum_id}{Vector of individual stratum IDs. nx1. Default is \code{NULL},
indicating no stratification.}

\item{V_data}{Dataframe of additional regression covariates. nxq. Factor
covariates must be converted to factors. If \code{NULL} (default), no additional
covariates are to be included. All variables in \code{glm_form} must
be found in \code{V_data}.}

\item{run_sampler}{String specifying which sampler(s) should be run. Must be
one of \code{"both"} (default), \code{"fixed"}, or \code{"adapt"}.}

\item{glm_form}{String specifying formula for probit regression, excluding
outcome and latent class. For example, \code{"~ 1"} for the model with only
latent class as covariates. All variables in \code{glm_form} must be found in \code{V_data}.
Do not specify interaction terms for latent class by additional covariates,
as these terms are already included.}

\item{K_max}{Upper limit for number of classes. Default is 30.}

\item{adapt_seed}{Numeric seed for adaptive sampler. Default is \code{NULL}.}

\item{class_cutoff}{Minimum class size proportion when determining number of
classes in adaptive sampler. Default is 0.05.}

\item{alpha_adapt}{Adaptive sampler hyperparameter for prior for class
membership probabilities \eqn{\pi}. Default is \code{NULL} and default values are
used. If specified, must be (\code{K_max})x1.}

\item{eta_adapt}{Adaptive sampler hyperparameter for prior for item consumption
level probabilities \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class
\eqn{k}, assumed to be the same across classes. Default is \code{NULL} and default
values are used. If specified, must be JxR, where J is the number of exposure
items and R is the maximum number of categories for the exposure.}

\item{alpha_fixed}{Fixed sampler hyperparameter for prior for \eqn{\pi}. Default is
\code{NULL} and default values are used. If specified, must be (\code{K_fixed})x1.}

\item{eta_fixed}{Fixed sampler hyperparameter for prior for item consumption
level probabilities \eqn{\theta_{jk\cdot}} for each item \eqn{j} and class
\eqn{k}, assumed to be the same across classes. Default is \code{NULL} and default
values are used. If specified, must be JxR, where J is the number of exposure
items and R is the maximum number of categories for the exposure.}

\item{K_fixed}{True number of classes, if known. Default is \code{NULL}, as it is
not necessary for the adaptive sampler. If bypassing the adaptive sampler and
running the fixed sampler directly, need to specify a value here.}

\item{fixed_seed}{Numeric seed for fixed sampler. Default is \code{NULL}.}

\item{n_runs}{Number of MCMC iterations. Default is 20000.}

\item{burn}{Number of MCMC iterations to drop as a burn-in period. Default is 10000.}

\item{thin}{Thinning factor for MCMC iterations. Default is 5.}

\item{ci_level}{Confidence interval level for probit regression coefficient
estimates. Default is \code{0.95}.}

\item{save_res}{Boolean specifying if results should be saved. Default = \code{TRUE}.}

\item{save_path}{String specifying directory and file name to save results. Default is \code{NULL}.}
}
\value{
If the fixed sampler is run, returns list \code{res} containing:
\describe{
\item{\code{estimates}}{List of posterior model results}
\item{\code{runtime}}{Total runtime for model}
\item{\code{data_vars}}{List of data variables used}
\item{\code{MCMC_out}}{List of full MCMC output}
\item{\code{post_MCMC_out}}{List of MCMC output after relabeling}
\item{\code{K_fixed}}{Number of classes used for the fixed sampler}
\item{\code{K_MCMC}}{If \code{K_fixed = NULL} and the adaptive sampler is run,
output list also contains MCMC output for the number of classes with size
greater than \code{class_cutoff} for each iteration}
}
If \code{save_res = TRUE} (default), also saves \code{res} as
\verb{[save_path]_wolca_results.RData}.

If only the adaptive sampler is run (i.e., \code{run_sampler} = \code{"adapt"}), returns
list \code{res} containing:
\describe{
\item{\code{MCMC_out}}{List of full MCMC output}
\item{\code{K_fixed}}{Number of classes used for the fixed sampler,
obtained from the adaptive sampler}
\item{\code{K_MCMC}}{Adaptive sampler MCMC output for K}
}
If \code{save_res = TRUE} (default), also saves \code{res} as
\verb{[save_path]_wolca_adapt.RData}.
}
\description{
\code{wolca} runs a two-step model with an unsupervised weighted overfitted latent
class analysis (WOLCA) in the first step and saves and returns the results.
}
\details{
\code{wolca} is a two-step approach that runs an unsupervised WOLCA in the first
step to derive latent class patterns and subsequently treats the class
assignments as fixed and includes them as covariates in a frequentist
survey-weighted probit regression model that uses an asymptotic sandwich
estimator for variance estimation.

By default, the function will run both samplers for the WOLCA step, running
the adaptive sampler first to determine the number of latent classes, and
then using the determined number of latent classes to run the fixed sampler
for parameter estimation of the unsupervised WOLCA, followed by the
survey-weighted regression model. If the number of latent classes is already
known and only the fixed sampler is to be run, specify \code{"fixed"} for the
\code{run_sampler} argument and specify a number for \code{K_fixed}. Id only the
adaptive sampler is to be run, specify \code{"adapt"} for the \code{run_sampler}
argument. Use \code{adapt_seed} (default is \code{NULL}) to specify a seed
for the adaptive sampler, and use \code{fixed_seed} (default is \code{NULL}) to specify
a separate seed for the fixed sampler.

\code{x_mat} is an nxJ matrix with each row corresponding to the J-dimensional
categorical exposure for an individual. If there is no clustering present,
\code{cluster_id} should be set to the individual IDs. \code{V_data} includes all
covariates to include in the probit regression other than latent class.
\code{K_max} is the maximum number of latent classes allowable, to
be used for the overfitted latent class model if the adaptive sampler is run.
\code{class_cutoff} is the minimum size of each class as a proportion of the
population, used when determining the number of latent classes.

To save results, set \code{save_res = TRUE} (default) and \code{save_path} to a string
that specifies both the location and the beginning of the file name
(e.g., "~/Documents/run"). The file name will have "_wolca_adapt.RData" or
"_wolca_results.RData" appended to it.

If hyperparameters for the adaptive or fixed sampler are left as \code{NULL}
(default), the following default values are used. Let \eqn{K} refer to
\code{K_max} for the adaptive sampler and \code{K_fixed} for the fixed sampler.
For \eqn{\pi}, a Dirichlet prior with hyperparameter \eqn{\alpha = 1/K} for
each component. For \eqn{\theta_{jk\cdot}}, a Dirichlet prior with
hyperparameter  \eqn{\eta_j} equal to \code{rep(1, R_j)} where \code{R_j} is the number
of categories for exposure item j. If \code{R_j < R}, the remaining categories have
hyperparameter set to 0.01. This is done independently for each exposure item j
and is assumed to be the same across latent classes. Note that hyperparameters
for the fixed sampler should probably only be specified if running the
fixed sampler directly, bypassing the adaptive sampler.
}
\examples{
# Load data and obtain relevant variables
data("sim_data")
data_vars <- sim_data
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
n <- dim(x_mat)[1]                   # Number of individuals

# Probit model only includes latent class
V_data <- NULL # Additional regression covariates
# Survey-weighted regression formula
glm_form <- "~ 1"

# Run wolca
res <- wolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt, 
       cluster_id = cluster_id, stratum_id = stratum_id, V_data = V_data, 
       run_sampler = "both", glm_form = glm_form, adapt_seed = 1, 
       n_runs = 50, burn = 25, thin = 1, save_res = FALSE)

}