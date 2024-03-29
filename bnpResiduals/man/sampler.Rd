% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bnpResiduals_functions.R
\name{sampler}
\alias{sampler}
\title{Sample all the parameters defining the model}
\usage{
sampler(
  X,
  Y,
  initializer,
  iter,
  burn_in = 10000,
  thin = 50,
  Verbose = T,
  log_rate = T,
  sigma = 0
)
}
\arguments{
\item{X}{design matrix}

\item{Y}{response variable}

\item{initializer}{list of initial values - can be obtained as the output of the function initializer()}

\item{iter}{number of iterations to consider}

\item{burn_in}{number of initial sampled valued to discard}

\item{thin}{thinning value}

\item{Verbose}{boolean value which allows to print information about the running of the algorithm}

\item{log_rate}{compute all the probabilities in log-scale (default: log_rate=T)}

\item{sigma}{parameter of the Pitman Yor process - default: sigma=0 (Dirichlet process case)}
}
\description{
This function allows you to obtain a Markov Chain with stationary distribution equal to the posterior of the parameters sampled, as explained in Ascolani, Ghidini (year)
}
\keyword{sampler}
