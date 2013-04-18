\name{salmodll}
\alias{salmodll}
\title{Call SALMO Shared Library}
\usage{
  salmodll(cfunc, nOfVar, cc, pp, uu, xx, pm, mx)
}
\arguments{
  \item{cfunc}{string, name of the C function to be called}

  \item{nOfVar}{vector with number of variables}

  \item{cc}{vector of constants for the model}

  \item{pp}{matrix of phytoplankton parameters}

  \item{uu}{input vector (environmental conditions)}

  \item{xx}{state vector}

  \item{pm}{parameters for macrophytes}

  \item{mx}{states of macrophytes (ReacTran format)}
}
\value{
  list, containing the derivatives as first element
}
\description{
  This function calls the shared library of SALMO (version
  with macrophyte coupling).
}

