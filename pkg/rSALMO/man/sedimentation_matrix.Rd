% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/sedimentation_matrix.R
\name{sedimentation_matrix}
\alias{sedimentation_matrix}
\title{Create Sedimentation and Migration Matrix for SALMO-1D}
\usage{
sedimentation_matrix(parms, nstates = NULL, nphy = NULL, depth = NULL,
  focussing = c(1, 1))
}
\arguments{
\item{parms}{SALMO parameter list}

\item{nstates}{total number of state variables  (deprecated !!!)}

\item{nphy}{number of phytoplankton groups   (deprecated !!!)}

\item{depth}{vector of individual depths (m)  (deprecated !!!)}

\item{focussing}{optional two-valued vector for a linear sediment focussing
heuristics (e.g. \code{c(1, 2)} - experimental!)}
}
\value{
matrix with sedimentation / migration velocities
}
\description{
Note: state variable indexes still hard coded
}
\details{
ToDo: pass parameters regularly to this function
}
\examples{
nstates <- 8 # nOfVar$numberOfStates
nphy    <- 3 # nOfVar$numberOfAlgae
parms   <- get_salmo_parms(nlayers=140, macrophytes=TRUE)
depth  <- seq(0, 70, 0.5) # unique(inputs[,3])
vmat   <- sedimentation_matrix(parms, nstates, nphy, depth)
}

