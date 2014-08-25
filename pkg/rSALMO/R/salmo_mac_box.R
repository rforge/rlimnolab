#' Derivatives of SALMO with macrophytes
#' 
#' Main model function to be called from ode solvers (1 box version with macrophytes).
#' 
#' 
#' @param time simulation time (days)
#' @param states state vector in correct order
#' @param parms list containing constant model parameters
#' @param inputs input vector (environmental conditions)
#' @return list, first element contains the derivatives, other elements can
#' contain optional outputs
#'
#' @rdname salmo_mac_box
#' @export salmo_mac_1box
salmo_mac_1box <- function(time, states, parms, inputs) {

  cat("time=", time, "\n")
  #print(states)
  if(any(states < 0)) {
    id <- which(states < 0)
    print(paste(time, " warning a state went negative and was set to 1e-7 ", id))
  }
  states <- ifelse(states < 0, 1e-7, states)

  with(parms, {
    # todo: make numbers flexible
    nstates <- 11
    nlayers <- 1

    y.salmo <- states[1:11]   # plankton submodel
    y.macro <- states[12:15]  # macrophyte submodel

    forc <- approxTime1(inputs$forcings, time)
    ## only required for systems > 0D
    #attr(forc, "colnames") <- attr(inputs$forcings, "colnames") # thpe: tricky :-(

    ## thpe: hard-coded indices "9" should be avoided. Use indirect variables
    #itemp <- which(attr(forcings, "colnames") == "temp")

    #temp       <- forc[itemp + (0:(nlayers - 1) * ni)] #layer temperature
    #zmixret    <- calczmix(temp, depths)
    #cc["Zres"] <- zmixret$zres          # set resuspension depth to mixing depth

    ## call SALMO core for every layer and returns derivatives for all layers
    dy.salmo <- numeric(length(y.salmo))

    ## reorder states for use with SALMO
    y.tmp  <- arrangeLayerWise(y.salmo, nOfVar["numberOfStates"],
                               nOfVar["numberOfLayers"])
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

    #aFunVegSed <- 1/c(1, salmo$aFunVegRes)
    aFunVegSed <- 1/salmo$aFunVegRes

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

    vsedi   <- inputs$vsedi   * aFunVegSed
    #vmatsedi   <- inputs$vmatsedi * aFunVegSed

    # thpe: transport for a single box
    # .... does not need ReacTran
    dtransport <- - vsedi[1,] * y.salmo / forc["dz"]

    ## this function would require matrices
    #dtransport <- transport(y.salmo, forc, parms, zmixret$idzmix, zmixret$zres,
    #  vmat, vmatsedi, time)

    ## state equation
    dy.salmo    <- dreaction[1:11] + dtransport + dmacro

    ## oxygen ==========
    #ynames <- names(y.salmo)
    #iO2 <-  which(ynames == "O")

    # the same "hard coded" for testing performance
    iO2 <- 8

    dy.salmo[iO2] <-  o2sat(forc["temp"]) - y.salmo[iO2]

    list(c(dy.salmo, macro$dy.macro))
  }) # end with
}
