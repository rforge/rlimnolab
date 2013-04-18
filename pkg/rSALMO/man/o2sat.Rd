\name{o2sat}
\alias{o2sat}
\title{Oxygen Saturation in Fresh Water}
\usage{
  o2sat(T)
}
\arguments{
  \item{T}{vector of temperatures (in deg C)}
}
\description{
  Simple formula for oxygen saturation (in mg/L) from
  temperature (in deg C).
}
\details{
  This is a simple function can also be used to test
  availability of the SALMO shared library
}
\examples{
o2sat(c(10, 20))
}
\seealso{
  \code{\link[marelac]{gas_O2sat}} in package marelac
}

