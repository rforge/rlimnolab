% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/salmo_1d.R
\name{salmo_1d}
\alias{salmo_1d}
\title{Derivatives of SALMO}
\usage{
salmo_1d(time, states, parms, inputs, ndx, forcingfun = NULL)
}
\arguments{
\item{time}{simulation time (days)}

\item{states}{state vector in correct order}

\item{parms}{list containing constant model parameters}

\item{inputs}{input vector (environmental conditions)}

\item{ndx}{hashtable (environment) of indexes and counters}

\item{forcingfun}{function that returns time dependent input data
  in appropriate order ........}
}
\value{
list, first element contains the derivatives, other elements can
contain optional outputs
}
\description{
Main model function to be called from ode solvers (1D version without macrophytes).
}

