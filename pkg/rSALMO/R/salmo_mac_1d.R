#' Derivatives of SALMO with macrophytes
#' 
#' Main model function to be called from ode solvers (1D version with macrophytes).
#' 
#' 
#' @param time simulation time (days)
#' @param states state vector in correct order
#' @param parms list containing constant model parameters
#' @param inputs input vector (environmental conditions)
#' @param ndx hashtable (environment) of indexes and counters 
#' @param forcingfun function that returns time dependent input data
#'   in appropriate order ........
#'   
#' @return list, first element contains the derivatives, other elements can
#' contain optional outputs
#'
#' @rdname salmo_mac_1d
#' @export salmo_mac_1d
salmo_mac_1d <- function(time, states, parms, inputs, ndx, forcingfun=NULL) {
  #cat("time = ", time, "\n")

  if(any(states < 0)) {
    id <- which(states < 0)
    #print(paste(time, " warning a state went negative and was set to 1e-7 ", id))
    cat(state_names[id], "\n")
  }
  states <- ifelse(states < 0, 1e-7, states)
  
  with(parms, {
    y.salmo <- states[1:(nlayers * nstates)]                   # plankton submodel
    y.macro <- states[(nlayers * nstates + 1):length(states)]  # macrophyte submodel
    
    ## forcings should be a function 'forcingfun'
    ## otherwise it has to be a full matrix
    if (is.null(forcingfun)) {
    forc <- approxTime1(inputs$forcings, time)
    attr(forc, "colnames") <- attr(inputs$forcings, "colnames") # thpe: tricky :-(
    } else {
      forc <- forcingfun(time)
    }
    
    itemp  <- ndx$itemp
    idepth <- ndx$idept
    temp       <- forc[itemp + (0:(nlayers - 1) * ni)] #layer temperature

    depth <- depths # thpe: fixme !!!
    zmixret    <- calczmix(temp, depths)
    # test test test
    ## optionally write calculated depths to log file
    if (syslog) cat(time, "\t", zmixret$idzmix, "\t", zmixret$zres, "\n", file="logfile.log", append = TRUE)
    
    ## limit resuspension depth by a given maximum value
    if (zmixret$zres > zresmax) {
      ## which.min works opposit to which here (finds 1st FALSE)
      idzmix <- which.min(depth <= zresmax) 
      zmixret <- list(idzmix=idzmix, zres = zresmax)
    }
    #cat(unlist(zmixret), "\n")

    cc["Zres"] <- zmixret$zres          # set resuspension depth to mixing depth
 
    ## call SALMO core for every layer and returns derivatives for all layers  
    dy.salmo <- numeric(length(y.salmo))
    ## reorder states for use with SALMO
    #y.tmp  <- arrangeLayerWise(y.salmo, nstates, nlayers)
    y.tmp <- as.vector(t(matrix(y.salmo, nrow=nlayers, ncol=nstates))) # layer wise
    
    salmo <- .C(
      "MReaktion",
      as.integer(nOfVar),
      as.double(cc),
      as.double(pp),
      forc       = as.double(forc),
      y.tmp      = as.double(y.tmp),    # contains modified forcings
      dy.salmo   = as.double(dy.salmo),
      as.double(pp_ma[,1]),             # only 1st macrophyte group implemented yet
      y.macro    = as.double(y.macro),
      aFunVegRes = as.double(numeric(nlayers))
    )
    
    aFunVegSed <- 1/c(1, salmo$aFunVegRes)
    
    # thpe: returned from from salmodll:
    #dreaction <- arrangeStateWise(salmo$dy.salmo, nstates, nlayers)
    dreaction <- as.vector(t(matrix(salmo$dy.salmo, nrow=nstates, ncol=nlayers))) #state wise
    
    ## call macrophyte model
    #salmodmx <- numeric(length(y.salmo)) # thpe: redundant
    #salmodmx <- numeric(length(salmo$xx))
    dy.macro <- numeric(length(y.macro))

    ## call C function
    macro <- .C(
      "MacrophyteReaction",
      as.double(nOfVar_ma),
      as.double(nOfVar),
      as.double(pp_ma[,1]),  # only 1st macrophyte group implemented yet
      as.double(cc_ma),
      as.double(salmo$forc),
      as.double(time),
      as.double(salmo$y.tmp),
      dy.salmo = as.double(dy.salmo),
      as.double(y.macro),
      dy.macro = as.double(dy.macro)
    )
     
    #dmacro <- arrangeStateWise(macro$dy.salmo, nstates, nlayers)
    dmacro <- as.vector(t(matrix(macro$dy.salmo, nrow=nstates, ncol=nlayers))) #state wise

    #vmat       <- inputs$vmat     * aFunVegSed
    vmatsedi   <- inputs$vmatsedi * aFunVegSed
    
    #dtransport <- transport(y.salmo, forc, parms, zmixret$idzmix, zmixret$zres, 
    #  vmat, vmatsedi, time)

    dtransport <- transport(y.salmo, forc, parms, ndx, zmixret$idzmix, zmixret$zres,
      vmatsedi)
    
    ## state equation
    dy.salmo    <- dreaction + dtransport + dmacro 
      
    list(c(dy.salmo, macro$dy.macro))
  }) # end with
}
