#' Derivatives of SALMO
#'
#' Main model function to be called from ode solvers
#' (Fully mixed one-box version).
#'
#' @param time   simulation time (days)
#' @param x      state vector in correct order
#' @param p      list containing constant model parameters
#' @param inputs  input vector (environmental conditions)
#'
#' @return list, first element contains the derivatives,
#'               other elements can contain optional outputs

SALMO.1box <- function(time, x, p, inputs) {
  #cat(time, "\n")

  ## interpolate data (this is the slowest part of the simulation)
  uu <- approxTime1(inputs, time)

  ## call the C model core
  ret <- call_salmodll("SalmoCore", p$nOfVar, p$cc, p$pp, uu, x)

  ## return source - sink terms
  dx <-c(ret[[2]] - ret[[3]])

  ## fix oxygen balance at surface to saturated value
  dx[8] <-  o2sat(uu["temp"]) - x[8]
  #list(dx, dO = unname(dx[8]))
  list(dx, iin = unname(uu["iin"]))
}
