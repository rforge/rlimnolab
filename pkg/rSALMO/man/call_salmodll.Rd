% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/call_salmodll.R
\name{call_salmodll}
\alias{call_salmodll}
\title{Call SALMO Shared Library}
\usage{
call_salmodll(cfunc, nOfVar, cc, pp, uu, xx)
}
\arguments{
\item{cfunc}{string, name of the C function to be called}

\item{nOfVar}{vector with number of variables}

\item{cc}{vector of constants for the model}

\item{pp}{matrix of phytoplankton parameters}

\item{uu}{input vector (environmental conditions)}

\item{xx}{state vector}
}
\value{
list, contains (1) past state and derivatives separated into (2)
source and (3) sink term and (4) modified inputs
}
\description{
Call the shared library of SALMO
}

