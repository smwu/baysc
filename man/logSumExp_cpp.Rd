% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{logSumExp_cpp}
\alias{logSumExp_cpp}
\title{Apply log-sum-exp trick}
\usage{
logSumExp_cpp(x)
}
\arguments{
\item{x}{Row vector of input values}
}
\value{
Result of computing the logarithm of the sum of exponentials of the
input values.
}
\description{
\code{logSumExp_cpp} computes the logarithm of the sum of exponentials, adapted from
https://github.com/helske/seqHMM/blob/main/src/logSumExp.cpp.
}
\keyword{internal}
