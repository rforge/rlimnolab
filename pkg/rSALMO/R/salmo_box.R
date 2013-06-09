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

salmo_1box <- function(time, x, p, inputs) {
  # save a few frequently used criteria to temporary variables
  #xnames <- names(x)
  #O2 <-  which(xnames == "O")
  
  # the same "hard coded" for testing performance
  O2 <- 8
  
  #cat(time, "\n")

  ## interpolate data (this is the slowest part of the simulation)
  uu <- approxTime1(inputs, time)

  ## call the C model core
  ret <- call_salmodll("SalmoCore", p$nOfVar, p$cc, p$pp, uu, x)

  ## return source - sink terms
  dx <-c(ret[[2]] - ret[[3]])

  ## fix oxygen balance at surface to saturated value
  
  dx[O2] <-  o2sat(uu["temp"]) - x[O2]
  
  list(dx, iin = unname(uu["iin"]))
}


salmo_2box <- function(time, x, p, inputs) {
  cat(time, "\n")
  noi <- nOfVar["numberOfInputs"]
  nos <- nOfVar["numberOfStates"]
  
  
  ## constants; should be made flexible
  undt <- 1  # index of eps, modified by phytoplankton; can be used by physics
  uvol <- 2  # index of volume of the layer
  uiin <- 8  # index of light at bottom of box  
  
  ## interpolate data (this is the slowest part of the simulation)
  uu <- approxTime1(inputs, time)

  uuE <- uu[1:noi]            # epilimnion  
  uuH <- uu[(noi+1):(2*noi)]  # hypolimnion
   
  ## call the C model core for the upper layer (epilimnion)
  retE <- call_salmodll("SalmoCore", p$nOfVar, p$cc, p$pp, uuE, x[1:nos])

  ## modify uu so that it describes the bottom layer (hypolimnion)
  # ... to be implemented
  
  ## light on top of hypo is ligt at bottom of epi
  uuH[uiin] <- retE[[4]][uiin]
  
  ## call the C model core for the bottom layer
  if (uuH[uvol] > 1e-14) {
    retH <- call_salmodll("SalmoCore", p$nOfVar, p$cc, p$pp, uuH, x[(nos+1):(2*nos)])
  } else {
    retH <- retE  # fully mixed conditions 
    # !!! this works only at beginning of the simulation !!!
    #     because states have also to be reset resp. balanced 
  }
  
  ## ========= Transport ===========
  ## apply partial mixing of the layers
  # ... to be implemented
  
  
  ## return source - sink terms
  dxE <-c(retE[[2]] - retE[[3]])
  dxH <-c(retH[[2]] - retH[[3]])

  ## Set oxygen balance at surface to saturated value
  dxE[8] <-  o2sat(uuE["temp"]) - x[8]
    
  list(c(dxE, dxH), uuE[uiin], uuH[uiin])
}
