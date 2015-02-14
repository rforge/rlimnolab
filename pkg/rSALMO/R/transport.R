#' Calulate Vertical Transport for a 1D Lake Model
#' 
#' This function calculates vertical transport for a 1D lake model.
#' It is normally not called by the user.
#' 
#' 
#' @param x vector of state variables
#' @param forcings matrix of environmental forcing data
#' @param parms list of model parameters
#' @param idzmix redundant ?
#' @param zres resuspension depth
#' @param vmatsedi matrix conmtaining the vertical sedimentation velocities for the state variables
#' @return vector of the derivatives of the model
#'
#' 
#' @export transport
#'
transport <- function(x, forcings, parms, ndx, idzmix, zres, vmatsedi) {

  dxx <- numeric(length(x))
  cnames <- attr(forcings, "colnames")
  
  isf   <- ndx$isf
  itemp <- ndx$itemp
  iaver <- ndx$iaver
  ivol  <- ndx$ivol
  idz   <- ndx$idz
  ieddy <- ndx$ieddy
  
  iO2 <- ndx$iO
  ni <- ndx$ninputs
  
  ## calculate area from volume and depth
  mforc <- matrix(forcings, nrow = ni)
  Area <- mforc[ivol, ] / mforc[idz,]
  Area <- c(Area, Area[length(Area)])
  
  with(parms, {

    # thpe: one of the matrices was obsolete
    vmat <- vmatsedi
    ## no flux.down of phytoplankton in last layer,
    ## because sedimentation to sediment is handled in SalmoCore itself
    vmat[nrow(vmat),] <- 0 #rep(0, nstates)

    Vol  <- forcings[ivol + (0:(nlayers - 1) * ni)] # layer volume
    dz   <- forcings[idz  + (0:(nlayers - 1) * ni)] # layer height

    D <- forcings[ieddy  + (1:(nlayers - 1) * ni)]
    D    <- c(0, D , 0) #no diffusion at upper and lower bounderies out of the system

    ## transport loop over all state variables
    for(i in 1:nstates){
      ## increase diffusion for zooplankton by a constant factor due to own mobility
      #if (i == iZ) {
      #  Dori <- D
      #  D    <- Dzoo
      #}
      ## positions of the i-th state variable for all layers (ReacTran format)
      id <- (i + (i-1)*(nlayers-1)) : ((i + (i-1)*(nlayers-1)) + nlayers - 1)
      tx <- x[id] # the i-th state variable for all layers
      ## flux of oxygen from atmosphere into first layer
      fluxup <- ifelse(i == iO2, K2 * (o2sat(forcings[itemp]) - tx[1]) / dz[1], 0)
      
      ## optional heuristics: 
      ##  - increase vsink for cyanos in autumn
      ##  - D Zoo + 1.5
      ## ...
      # vmat[,6]          <- -0.15 #VMIG 
      
      ## transport of the i-th state variable
      dxx[id] <- tran.1D(tx, C.up = 0, C.down = 0, D = D, v = vmat[,i],
         flux.up = fluxup, A = Area, dx = dz)$dC
      #if (i == 6) D <- Dori
    }

    ## sedimentation into sediment (only calculated when external sedimentation
    ## is switched on / internal sedimentation is switched off)
    if(cc["SF"] == 0){
      depthsvec <- rep(depths, nstates)
      sedicheck <- ifelse(depthsvec > zres, TRUE, FALSE)
      #cat("zres=", zres, "\n")
      sedi <- as.vector(vmatsedi[2:nrow(vmatsedi),]) / rep(dz, nstates) * # vsink / dz
        rep(forcings[isf + (0:(nlayers-1)) * ni], nstates) *              # SF ## thpe: instead of vmatsedi ???!!!
        #(1 - rep(forcings[iaver + (0:(nlayers-1)) * ni], nstates)) *      # 1 - aver
        x * sedicheck                                                     # depth > zres

      ## thpe: add check for negative sedimentation
      if (any(sedi < 0)) {
        cat("warning: negative sedimentation!\n")
        sedi <- pmax(sedi, 0)
      }
      dxx <- dxx - sedi
    }
    ## return derivatives
    dxx
  })
}
