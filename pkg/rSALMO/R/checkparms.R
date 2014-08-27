#' Check SALMO Model Parameters
#' 
#' Consistency checks of model parameters and creation of internal counters.
#' (Consistency checks are not yet implemented)
#' 
#' @param parms SALMO parameter list
#' @param optpar list of additional 'global parameters' like depths or K2
#'   (shall be removed in the future)
#' @return the updated \code{parms} list
#'
#' @export check_salmo_parms
#'
check_salmo_parms <- function(parms, optpar) {
  ## do checks ...
  # - order of names
  # - length of vectors and matrices
  # - ...
  
  ## add several "global" parameters that frequently used
  parms$nstates  <- parms$nOfVar["numberOfStates"]
  parms$nlayers  <- parms$nOfVar["numberOfLayers"]
  parms$nphy     <- parms$nOfVar["numberOfAlgae"]
  parms$ni       <- parms$nOfVar["numberOfInputs"]
  parms$nspec    <- parms$nOfVar["numberOfStates"] + parms$nOfVar_ma["numberOfStates"] 
  parms <- c(parms, optpar)
  parms
}
