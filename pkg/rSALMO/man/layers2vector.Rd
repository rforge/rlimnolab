\name{arrangeStateWise}
\alias{arrangeLayerWise}
\alias{arrangeStateWise}
\alias{layers2vector}
\title{Convert Layer Structure to Vector}
\usage{
  arrangeStateWise(x, nstates, nlayers)

  arrangeLayerWise(x, nstates, nlayers)

  layers2vector(x)
}
\arguments{
  \item{x}{state vector resp. state matrix}

  \item{nstates}{number of state variables}

  \item{nlayers}{number of layers}
}
\value{
  rearranged vector of state variables
}
\description{
  These functions are used to reorganize the state vector
  state-wise resp. layer-wise.
}
\details{
  Functions \code{arrangeStateWise} and
  \code{arrangeLayerWise} rearrange state vectors from
  layer wise (ABCABCABC) to state wise (AAABBBCCC) and vice
  versa, function \code{layers2vector} creates a state
  vector from a matrix (layers x variables).
}

