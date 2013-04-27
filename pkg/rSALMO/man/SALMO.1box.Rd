\name{SALMO.1box}
\alias{SALMO.1box}
\title{Derivatives of SALMO}
\usage{
  SALMO.1box(time, x, p, inputs)
}
\arguments{
  \item{time}{simulation time (days)}

  \item{x}{state vector in correct order}

  \item{p}{list containing constant model parameters}

  \item{inputs}{input vector (environmental conditions)}
}
\value{
  list, first element contains the derivatives, other
  elements can contain optional outputs
}
\description{
  Main model function to be called from ode solvers (Fully
  mixed one-box version).
}

