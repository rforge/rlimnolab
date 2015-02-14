#' Initialize Index Values of Parameter List
#' 
#' Add precalculated index values to the parameter list, to avoid repeated
#' searching.
#' 
#' 
#' @param parms list containing constant model parameters
#' 
#' @return environment (hash table) that indexes for inputs and states
#'
#' @rdname init_salmo_integers
#' @export init_salmo_integers
#'
init_salmo_integers <- function(parms) {

  ndx <- list()
  
  ## indexes for state variables
  states <- salmo_state_names(nlayers = 1)
  ssn <- data.frame(name = states, ndx = 1:length(states))
  ssn$iName <- paste("i", ssn$name, sep="")
  ndx[ssn$iName] <- ssn$ndx
  
  ## indexes for state variables
  inputs <- salmo_input_names()
  ssn <- data.frame(name = inputs, ndx = 1:length(inputs))
  ssn$iName <- paste("i", ssn$name, sep="")
  ndx[ssn$iName] <- ssn$ndx
  
  ## add counters
  ndx$ninputs  <- length(inputs)
  ndx$nstates  <- parms$nOfVar["numberOfStates"]
  ndx$nlayers  <- parms$nOfVar["numberOfLayers"]
  ndx$nphy     <- parms$nOfVar["numberOfAlgae"]
  ndx$ni       <- parms$nOfVar["numberOfInputs"]
  
  ## use an environment as hash table
  as.environment(ndx)
}