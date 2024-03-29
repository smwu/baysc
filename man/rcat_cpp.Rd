% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rcat_cpp}
\alias{rcat_cpp}
\title{Draw from a Categorical distribution}
\usage{
rcat_cpp(probs)
}
\arguments{
\item{probs}{Row vector of category event probabilities. Length determines
the number of categories for the Categorical distribution.}
}
\value{
Integer specifying the category (1-based indexing) of a single draw
from the specified Categorical distribution
}
\description{
Draw from a Categorical distribution
}
\keyword{internal}
