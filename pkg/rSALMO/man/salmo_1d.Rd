\name{salmo_1d}
\alias{salmo_1d}
\title{Derivatives of SALMO}
\usage{
  salmo_1d(time, states, parms, inputs, forcingfun = NULL)
}
\arguments{
  \item{time}{simulation time (days)}

  \item{states}{state vector in correct order}

  \item{parms}{list containing constant model parameters}

  \item{inputs}{input vector (environmental conditions)}
}
\value{
  list, first element contains the derivatives, other
  elements can contain optional outputs
}
\description{
  Main model function to be called from ode solvers (1D
  version without macrophytes).
}

