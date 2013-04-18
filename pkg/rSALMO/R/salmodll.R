#' Call SALMO Shared Library
#' 
#' This function calls the shared library of SALMO 
#' (version with macrophyte coupling).
#' 
#' 
#' @param cfunc  string, name of the C function to be called
#' @param nOfVar vector with number of variables
#' @param cc     vector of constants for the model
#' @param pp     matrix of phytoplankton parameters
#' @param uu     input vector (environmental conditions)
#' @param xx     state vector
#' @param pm     parameters for macrophytes
#' @param mx     states of macrophytes (ReacTran format)
#' 
#' @return list, containing the derivatives as first element
salmodll <-
function(
  cfunc,         # name of the C function to be called
  nOfVar,        # number of variables
  cc,            # vector of constants for the model
  pp,            # matrix for phytoplankton parameters
  uu,            # vector of boundary conditions
  xx,            # state variables
  pm,            # parameters for macrophytes
  mx             # states of macrophytes (ReacTran format)
){
  ## create empty data structure for derivatives
  dxx <- numeric(length(xx))
  ## sort order of states for use with SALMO
  xx  <- sortSalmo(xx, nOfVar["numberOfStates"], nOfVar["numberOfLayers"])
  ## call C function
  ## thpe -> rene: variablen brauchen nur dann einen namen, wenn man sie 
  ##               wieder zurueckgeben will. 
  ret <- .C(
    cfunc,
    as.integer(nOfVar),
    cc = as.double(cc),
    pp = as.double(pp),
    uu = as.double(uu),
    xx = as.double(xx),
    dxx = as.double(dxx),
    pm = as.double(pm),
    mx = as.double(mx),
    aFunVegRes = as.double(numeric(nOfVar["numberOfLayers"]))
  )
  ## sort derivatives for use with ReacTran, 
  ## Attention: inputs uu and states xx are still in SALMO format 
  ##            for use with macrophyte module
  list(arrangeStateWise(ret$dxx, nOfVar["numberOfStates"], 
    nOfVar["numberOfLayers"]), ret$uu, ret$xx, ret$aFunVegRes)
  ## SALMO format of derivatives
  #  list(ret$dxx, ret$uu, ret$xx)
}
