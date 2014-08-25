#' Derivatives of SALMO with macrophytes
#' 
#' Main model function to be called from ode solvers (1D version with macrophytes).
#' 
#' 
#' @param time simulation time (days)
#' @param states state vector in correct order
#' @param parms list containing constant model parameters
#' @param inputs input vector (environmental conditions)
#' @return list, first element contains the derivatives, other elements can
#' contain optional outputs
#'
#' @rdname salmo_mac_1d
#' @export salmo_mac_1d
salmo_mac_1d <- function(time, states, parms, inputs) {
  cat("time = ", time, "\n")

  if(any(states < 0)) {
    id <- which(states < 0)
    #print(paste(time, " warning a state went negative and was set to 1e-7 ", id))
    cat(state_names[id], "\n")
  }
  states <- ifelse(states < 0, 1e-7, states)
  
  with(parms, {
    y.salmo <- states[1:(nlayers * nstates)]                   # plankton submodel
    y.macro <- states[(nlayers * nstates + 1):length(states)]  # macrophyte submodel
    
    forc <- approxTime1(inputs$forcings, time)
    attr(forc, "colnames") <- attr(inputs$forcings, "colnames") # thpe: tricky :-(
    
    #cat(attr(forc, "colnames"), "\n\n")
    
    ## thpe: hard-coded indices "9" should be avoided. Use indirect variables
    itemp <- which(attr(forcings, "colnames") == "temp")
    #temp      <- forc[9 + (0:(nlayers - 1) * ni)] #layer temperature 
    temp       <- forc[itemp + (0:(nlayers - 1) * ni)] #layer temperature 
    zmixret    <- calczmix(temp, depths)
    cc["Zres"] <- zmixret$zres          # set resuspension depth to mixing depth
    
    ## todo: 
    ##   make order of function arguments more logical and more consistent
    ##   - e.g. scalars first (e.g. t, ...
    ##   - or: like deSolve (t, x, p)
    ##              
  
    ## call SALMO core for every layer and returns derivatives for all layers  
    dy.salmo <- numeric(length(y.salmo))
    ## reorder states for use with SALMO
    y.tmp  <- arrangeLayerWise(y.salmo, nOfVar["numberOfStates"], nOfVar["numberOfLayers"])
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
    dreaction <- arrangeStateWise(salmo$dy.salmo, nstates, nlayers)
    
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
     
    dmacro <- arrangeStateWise(macro$dy.salmo, nstates, nlayers)

    vmat       <- inputs$vmat     * aFunVegSed
    vmatsedi   <- inputs$vmatsedi * aFunVegSed
    
    dtransport <- transport(y.salmo, forc, parms, zmixret$idzmix, zmixret$zres, 
      vmat, vmatsedi, time)
    
    ## state equation
    dy.salmo    <- dreaction + dtransport + dmacro 
      
    list(c(dy.salmo, macro$dy.macro))
  }) # end with
}
