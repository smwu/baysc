% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_functions.R
\name{get_param_mcmc}
\alias{get_param_mcmc}
\title{Obtain parameter MCMC iterations}
\usage{
get_param_mcmc(res)
}
\arguments{
\item{res}{An object of class \code{"swolca"} or \code{"wolca"}, resulting from a call
to \code{\link[=swolca]{swolca()}}, \code{\link[=swolca_var_adjust]{swolca_var_adjust()}}, \code{\link[=wolca]{wolca()}}, or \code{\link[=wolca_var_adjust]{wolca_var_adjust()}}.}
}
\value{
Returns list \code{param_mcmc} containing the following:
\describe{
\item{\code{pi_mcmc}}{MxK dataframe of the class membership probability
parameters, \eqn{\pi}, where M is the number of iterations and K is the
number of latent classes.}
\item{\code{theta_mcmc}}{Mx(JxKxR) dataframe of the item level probability
parameters, \eqn{\theta}, where J is the number of items and R is the
maximum number of item levels.}
}
If \code{res} is a \code{"swolca"} object, \code{param_mcmc} also contains \code{xi_mcmc}, a
Mx(KxQ) dataframe of the regression parameters, \eqn{\xi}, where Q is the
number of covariates, excluding latent class indicators, in the regression.
}
\description{
\code{get_param_mcmc} extracts the MCMC iteration values for all parameters
obtained through the MCMC sampler.
}
\examples{
data(run_nhanes_swolca_results)
param_mcmc <- get_param_mcmc(res = run_nhanes_swolca_results)

}
\seealso{
\code{\link[=summarize_res]{summarize_res()}}
}
