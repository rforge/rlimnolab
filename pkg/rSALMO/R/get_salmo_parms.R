#' Get SALMO Parameters
#'
#' Read default model parameters and construct parameter list.
#'
#'
#' @param nlayers integer specifying number of layers in 1D case,
#'   must be 1 even in the two-box setting.
#' @param macrophytes TRUE if macrophyte parametes shoudl be included
#'
#' @return a list with all required model parameters
#'
#' @export get_salmo_parms


get_salmo_parms <- function(nlayers = 1, macrophytes = FALSE) {

  ## unique SALMO parameters, vector
  data(parms_salmo_other)
  cc        <- parms_salmo_other$value
  names(cc) <- parms_salmo_other$id

  ## phytoplankton parameters, matrix with 1 column per phytopl. species
  data(parms_salmo_phyto)
  pp        <- parms_salmo_phyto[,-1]
  row.names(pp) <- parms_salmo_phyto[,1]
  pp <- as.matrix(pp)



  ## NOFVAR
  nOfVar <- c(
    numberOfInputs      = 22,  # Check: 21 for box versions?
    numberOfOutputs     = 14,
    numberOfStates      = 11,
    numberOfParameters  = NA,
    numberOfAlgae       = 3,
    numberOfLayers      = as.integer(nlayers),
    numberOfTributaries = 1,
    numberOfOutlets     = 1,
    timestep            = 1   # check: essential / obsolete / for compatibility only ?
  )
  
  #if (nlayers > 1) nOfVar$numerOfInputs <- 22 # contains SF column

  ## important!
  nOfVar["numberOfParameters"] <- length(pp) / nOfVar["numberOfAlgae"]

  ret <- list(
    pp = pp,
    cc = cc,
    nOfVar = nOfVar
  )

  ## todo: macrophyte parameters ...
  if (macrophytes) {
    ## unique macrophytes parameters, vector
    data(parms_mac_other)
    cc_ma        <- parms_mac_other$value
    names(cc_ma) <- parms_mac_other$id

    ## phytoplankton parameters, matrix with 1 column per phytopl. species
    data(parms_mac_plants)
    pp_ma         <- parms_mac_plants[,-1]
    row.names(pp_ma) <- parms_mac_plants[,1]
    pp_ma <- as.matrix(pp_ma)

    data(parms_mac_ctrl)
    nOfVar_ma        <- parms_mac_ctrl$value
    names(nOfVar_ma) <- parms_mac_ctrl$id
    ## total number of layers
    nOfVar_ma["numberOfLayers"] <- nOfVar["numberOfLayers"]
    ret <- list(
      pp = pp,
      cc = cc,
      nOfVar = nOfVar,
      pp_ma = pp_ma,
      cc_ma = cc_ma,
      nOfVar_ma = nOfVar_ma
    )
  }

  return(ret)
}
