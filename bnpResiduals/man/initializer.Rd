% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bnpResiduals_functions.R
\name{initializer}
\alias{initializer}
\title{Initialize the blocked Gibbs sampler algorithm}
\usage{
initializer(Z, p, comp, N, mu0, sigma0, complete_matrix, b0 = NULL, B0 = NULL)
}
\arguments{
\item{Z}{matrix containing the unique values of the triplets (Z1, Z2, Z3), each one associated with a label in 1:N}

\item{p}{distribution of the N labels}

\item{comp}{vector of dimension n (= number of observations), with possible entries +1 or -1, denoting the element of the mixture from which the residual has been sampled from}

\item{N}{number of elements in the truncation of the Pitman-Yor process}

\item{mu0}{hyperparameter for the baseline measure (mean of the Truncated Normal)}

\item{sigma0}{hyperparameter for the baseline measure (variance of the Truncated Normal)}

\item{complete_matrix}{matrix of dimension nx3 containing all the triplets associated to each observation}

\item{b0}{mean vector for the gaussian prior of beta - to specify just in case of linear models}

\item{B0}{covariance matrix for the gaussian prior of beta - to specify just in case of linear models}
}
\description{
This function allows you to define the initialization values for the blocked Gibbs sampler algorithm, with the quantities illustrated in Ascolani, Ghidini (2021).
All the quantities in input can be random, since they should not impact the final linear model estimation.
}
\examples{
initializer <- initializer(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix)
}
\keyword{initializer}
