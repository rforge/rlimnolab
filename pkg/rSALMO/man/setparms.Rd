\name{setparms}
\alias{setparms}
\title{Set Model Parameters}
\usage{
  setparms(obj, pnames, values, col = 1, slot = NA)
}
\arguments{
  \item{obj}{data object containing model parameters
  (currently only vector and matrix implemented)}

  \item{pnames}{vector of parameter names that shoudl be
  replaced}

  \item{values}{vector of replacement values}

  \item{col}{column in matrix-like parameter objects}
}
\value{
  the manipulated \code{obj} (vector, matrix or list)
}
\description{
  Robustly replace model parameters by ignoring
  capitalisation and order of replacment.
}
\examples{
data(cc)
cc <- setparms(cc, c("YZP","dummy", "EPSMIN", "dtt"), c(0.7, 66, 0.77, 0.1))
cc

data(pp)
pp <- setparms(pp, c("TOPTX", "vs"), c(28, 0.22), c(1, 3))
# pp <- setparms(pp, "KP", c(1.4, 1.5, 1.6), 1:3) # not yet possible
}

