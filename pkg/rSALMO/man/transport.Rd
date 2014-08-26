\name{transport}
\alias{transport}
\title{Calulate Vertical Transport for a 1D Lake Model}
\usage{
  transport(x, forcings, parms, idzmix, zres, vmatsedi)
}
\arguments{
  \item{x}{vector of state variables}

  \item{forcings}{matrix of environmental forcing data}

  \item{parms}{list of model parameters}

  \item{idzmix}{redundant ?}

  \item{zres}{resuspension depth}

  \item{vmatsedi}{matrix conmtaining the vertical
  sedimentation velocities for the state variables}
}
\value{
  vector of the derivatives of the model
}
\description{
  This function calculates vertical transport for a 1D lake
  model. It is normally not called by the user.
}

