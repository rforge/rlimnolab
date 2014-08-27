\name{check_salmo_parms}
\alias{check_salmo_parms}
\title{Check SALMO Model Parameters}
\usage{
  check_salmo_parms(parms, optpar)
}
\arguments{
  \item{parms}{SALMO parameter list}

  \item{optpar}{list of additional 'global parameters' like
  depths or K2 (shall be removed in the future)}
}
\value{
  the updated \code{parms} list
}
\description{
  Consistency checks of model parameters and creation of
  internal counters. (Consistency checks are not yet
  implemented)
}

