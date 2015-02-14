#' Derivatives of SALMO
#' 
#' Main model function to be called from ode solvers (Fully mixed one-box
#' version).
#' 
#' 
#' @param time simulation time (days)
#' @param x state vector in correct order
#' @param p list containing constant model parameters
#' @param inputs input vector (environmental conditions)
#' @param ndx hashtable (environment) of indexes and counters 
#' @return list, first element contains the derivatives, other elements can
#' contain optional outputs
#'
#' @rdname salmo_box
#' @export salmo_1box

salmo_1box <- function(time, x, p, inputs, ndx) {
  #cat(time, "\n")

  ## interpolate data (this is the slowest part of the simulation)
  uu <- approxTime1(inputs, time)

  ## call the C model core
  ret <- call_salmodll("SalmoCore", p$nOfVar, p$cc, p$pp, uu, x)

  ## return source - sink terms
  dx <-c(ret[[2]] - ret[[3]])

  ## fix oxygen balance at surface to saturated value
  dx[ndx$iO] <-  o2sat(uu["temp"]) - x[ndx$iO]
    
  ## sedimentation
  # to be checked
  v     <- unname(c(0, 0, p$pp["VS", ], 0, p$cc["VD"], 0, 0, 0, 0))
  sed <- x * v / uu["dz"]
  
  list(dx - sed, iin = unname(uu["iin"]))
}


#' @rdname salmo_box
#' @export salmo_2box
salmo_2box <- function(time, x, p, inputs, ndx) {
  #cat("Salmo 2box time=", time, "\n")
  
  
  if(any(x < 0)) {
    id <- which(x < 0)
    #print(paste(time, " warning a state went negative and was set to 1e-7 ", id))
    #cat(x, "\n")
  }
  x <- ifelse(x < 0, 1e-7, x)
  
  noi <- ndx$ninputs
  nos <- ndx$nstates
  
  
  ## constants; should be made flexible
  #undt <- ndx$itime # 1  # index of eps, modified by phytoplankton; can be used by physics
  uvol <- ndx$ivol  # 2  # index of volume of the layer
  uiin <- ndx$iiin  # 8  # index of light at bottom of box  

  ## interpolate data (this is the slowest part of the simulation)
  uu <- approxTime1(inputs, time)
  
  uuE <- uu[1:noi]            # epilimnion  
  uuH <- uu[(noi+1):(2*noi)]  # hypolimnion

  
  stratified <- uuH[uvol] > 1e-14
    
  xE <- x[1:nos]
  xH <- x[(nos+1):(2*nos)]
   
  ## call the C model core for the upper layer (epilimnion)
  retE <- call_salmodll("SalmoCore", p$nOfVar, p$cc, p$pp, uuE, xE)


  ## light on top of hypo is ligt at bottom of epi
  uuH[uiin] <- retE[[4]][uiin]
  
  
  ## call the C model core for the bottom layer
  if (uuH[uvol] > 1e-14) {
    #cat("light", uuE[uiin], uuH[uiin], "\n")
    retH <- call_salmodll("SalmoCore", p$nOfVar, p$cc, p$pp, uuH, xH)
  } else {
    retH <- retE  # fully mixed conditions 
    # !!! this works only at beginning of the simulation !!!
    #     because states have also to be reset resp. balanced 
  }
  
  ## ========= Transport ===========
  ## apply partial mixing of the layers
  # ... to be implemented
  
  if (stratified) {
    # !! check this VERY carefully !!!
    trE <- (xH-xE) * uuE["ad"] / uuE["vol"]
    trH <- (xE-xH) * uuE["au"] / uuH["vol"]
    
    #cat("trE ", trE, "\n")
    #cat("trH ", trH, "\n")
  } else {
    trE <- trH <- 0
  }
  
  ## sedimentation
  #v     <- unname(c(0, 0, p$pp["VS", ], 0, p$cc["VD"], 0, 0, 0, 0))
  v <- rep(0, nos)
  ## sedimentation velocity VS of X1, X2 and D
  v[c(ndx$iX1, ndx$iX2, ndx$iX3, ndx$iD)] <- c(p$pp["VS",], p$cc["VD"])
  
  
  sedE <- xE * v / uuE["dz"]
  
  ## to follow the original: v * 2 or v * SF ?
  if (stratified) {
    sedH <- xH * v * p$cc["SF"] / uuH["dz"] 
  } else {
    sedH <- 0
    
  }
  
  
  # todo: zmig 
  
  ##      source - sink terms + transp - sed.loss + sed.from.epi * area
  dxE <-c(retE[[2]] - retE[[3]]) + trE - sedE 
  dxH <-c(retH[[2]] - retH[[3]]) + trH - sedH +  sedE *  uuE["aver"]

  ## Set oxygen balance at surface to saturated value
  dxE[ndx$iO] <-  o2sat(uuE["temp"]) - x[ndx$iO]
  
  ## correct also for "dummy hypo" if no stratification
  if (! stratified) {
    dxH[ndx$iO] <-  o2sat(uuE["temp"]) - x[ndx$iO + ndx$nstates] # ndx = 8, 19
  }
    
  list(c(dxE, dxH))#, uuE[uiin], uuH[uiin])
}
