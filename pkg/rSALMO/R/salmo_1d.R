#' Derivatives of SALMO
#' 
#' Main model function to be called from ode solvers (1D version without macrophytes).
#' 
#' 
#' @param time simulation time (days)
#' @param states state vector in correct order
#' @param parms list containing constant model parameters
#' @param inputs input vector (environmental conditions)
#' @param forcingfun function that returns time dependent input data
#'   in appropriate order ........
#'   
#' @return list, first element contains the derivatives, other elements can
#' contain optional outputs
#'
#' @export salmo_1d

salmo_1d <- function(time, states, parms, inputs, forcingfun=NULL) {

  if(any(states < 0)) {
    id <- which(states < 0)
    #print(paste(time, " warning a state went negative and was set to 1e-7 ", id))
    cat(state_names[id], "\n")
  }
  states <- ifelse(states < 0, 1e-7, states)

  with(parms, {
    y.salmo <- states[1:(nlayers * nstates)]                   # plankton submodel
    #y.macro <- states[(nlayers * nstates + 1):length(states)]  # macrophyte submodel

    ## forcings should be a function 'forcingfun'
    ## otherwise it has to be a full matrix
    if (is.null(forcingfun)) {
      cat("time = ", time, "\n")
      forc <- approxTime1(inputs$forcings, time)
      attr(forc, "colnames") <- attr(inputs$forcings, "colnames") # thpe: tricky :-(
    } else {
      forc <- forcingfun(time)
    }

    itemp  <- which(attr(forc, "colnames") == "temp")
    idepth <- which(attr(forc, "colnames") == "depth")
    temp       <- forc[itemp + (0:(nlayers - 1) * ni)] #layer temperature
    #depth      <- forc[idepth + (0:(nlayers - 1) * ni)] #layer depth
    depth <- depths # thpe: fixme !!!
    zmixret    <- calczmix(temp, depth)

    ## optionally write calculated depths to log file
    if (syslog) cat(time, "\t", zmixret$idzmix, "\t", zmixret$zres, "\n", file="logfile.log", append = TRUE)
    
    ## limit resuspension depth by a given maximum value
    if (zmixret$zres > zresmax) zmixret <- list(idzmix=41, zres = zresmax)
    

    cc["Zres"] <- zmixret$zres          # set resuspension depth to mixing depth

    ## call SALMO core for every layer and returns derivatives for all layers
    dy.salmo <- numeric(length(y.salmo))

    ## reorder states for use with SALMO
    y.tmp  <- arrangeLayerWise(y.salmo, nstates, nlayers)
    salmo <- .C(
      "RReaktion",
      as.integer(nOfVar),
      as.double(cc),
      as.double(pp),
      forc       = as.double(forc),
      y.tmp      = as.double(y.tmp),    # Note: contains modified forcings
      dy.salmo   = as.double(dy.salmo)
    )

    # thpe: returned from from salmodll:
    dreaction <- arrangeStateWise(salmo$dy.salmo, nstates, nlayers)

    ## dummy for macrophyte model (zero derivatives)
    #dy.macro <- numeric(length(y.macro))

    dtransport <- transport(y.salmo, forc, parms, zmixret$idzmix, zmixret$zres,
      inputs$vmatsedi)

    ## state equation
    dy.salmo    <- dreaction + dtransport #+ dmacro

    #list(c(dy.salmo, dy.macro))
    list(dy.salmo)
  }) # end with
}
