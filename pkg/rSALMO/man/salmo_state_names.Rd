\name{salmo_input_names}
\alias{salmo_input_names}
\alias{salmo_state_names}
\title{Column Names of SALMO State and Input Variables}
\usage{
  salmo_input_names()

  salmo_state_names(nlayers, macrophytes = FALSE)
}
\arguments{
  \item{macrophytes}{logical if state variable names for
  macrophytes should be included}

  \item{nlayers}{numeric, number of layers fro replicating
  the names}

  \item{parms}{list containing constant model parameters}
}
\value{
  vector of names for states resp inputs. contain optional
  outputs
}
\description{
  Column Names of SALMO State and Input Variables
}

