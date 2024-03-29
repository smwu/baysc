% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_data.R
\name{get_betas_c}
\alias{get_betas_c}
\title{Obtain betas matrix for generating categorical latent class assignment}
\usage{
get_betas_c(pi_mat, formula_c, V_unique)
}
\arguments{
\item{pi_mat}{Matrix where each row is the class assignment probabilities for
a level of the categorical covariate. HxK, where H is the number of levels of
the categorical covariate and K is the number of latent classes. Rows of
\code{pi_mat} must sum to 1.}

\item{formula_c}{String specifying formula for multinomial logistic
regression to create category latent class assignment C.}

\item{V_unique}{Dataframe with one column containing the unique values of
the categorical covariate specified in \code{formula_c}. If \code{formula_c = "~1"},
set \code{V_unique = NULL}.}
}
\value{
Returns \code{beta_mat} matrix of betas to be used in a multinomial logistic
regression to generate a categorical variable C. \code{beta_mat} has K rows and
number of columns equal to the number of levels in the categorical covariate.
}
\description{
Obtain matrix of betas that can be used to generate the categorical latent
class assignment variable C using a multinomial logistic regression where C
may depend on a categorical covariate such as the stratum variable S.
}
\examples{
## Example 1: latent class C depends on stratum variable S
# Number of latent classes and number of levels of S
K <- 3; H <- 2
# Formula specifying that C depends on S
formula_c <- "~ s_all"
# Dataframe with unique values of S
V_unique <- data.frame(s_all = as.factor(1:H))
# Matrix of class assignment probabilities for each level of S
pi_mat <- matrix(c(0.3, 0.5, 0.2,   # class membership probs for S=1
                   0.1, 0.6, 0.3),  # class membership probs for S=2
                 byrow = TRUE, nrow = H, ncol = K)
# Get matrix of betas for generating C
beta_mat_c <- get_betas_c(pi_mat = pi_mat, formula_c = formula_c, 
                          V_unique = V_unique)
beta_mat_c

## Example 2: latent class is generated independently of other variables
# Matrix of class assignment probabilities
pi_mat <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
formula_c <- "~ 1"
V_unique <- NULL
get_betas_c(pi_mat = pi_mat, formula_c = formula_c, V_unique = V_unique)

}
