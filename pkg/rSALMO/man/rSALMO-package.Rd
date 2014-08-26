\docType{package}
\name{rSALMO-package}
\alias{rSALMO}
\alias{rSALMO-package}
\title{Ecological Simulation of Lakes}
\description{
  \tabular{ll}{ Package: \tab rSALMO\cr Type: \tab
  Package\cr Version: \tab 0.1\cr Date: \tab 2013-04-27\cr
  License: \tab GPL (>= 2)\cr LazyLoad: \tab yes\cr }
}
\details{
  An R implementation of the SALMO lake model\cr (SALMO:
  Simulation by means of An analytical Lake MOdel)
}
\examples{
demo("salmo_1box")
demo("salmo_2box")
}
\author{
  Rene Sachse and Thomas Petzoldt (package)\cr

  Susanne Rolinski (C code of SALMO)\cr

  Juergen Benndorf + (system of equations of SALMO)\cr

  Maintainer: Thomas Petzoldt
  <thomas.petzoldt@tu-dresden.de>\cr
}
\references{
  see
  \url{http://hhbio.wasser.tu-dresden.de/projects/salmo/}
}
\seealso{
  \code{\link[deSolve:deSolve-package]{deSolve}}
}
\keyword{package}

