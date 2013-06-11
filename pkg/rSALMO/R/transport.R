transport <- function(x, forcings, parms, idzmix, zres, vmat, vmatsedi, time) {
  dxx <- numeric(length(x))
  cnames <- attr(forcings, "colnames")
  
  isf   <- which(cnames == "sf")
  itemp <- which(cnames == "temp")
  iaver <- which(cnames == "aver")
  ivol  <- which(cnames == "vol")
  idz   <- which(cnames == "dz")
  
  iO2 <- 8
  ni <- parms$nOfVar["numberOfInputs"]
  
  ## calculate area from volume and depth
  mforc <- matrix(forcings, nrow = ni)
  Area <- mforc[ivol, ] / mforc[idz,]
  Area <- c(Area, Area[length(Area)])

  with(parms, {

    Vol  <- forcings[ivol + (0:(nlayers - 1) * ni)] # layer volume
    dz   <- forcings[idz  + (0:(nlayers - 1) * ni)] # layer height

    ## diffusion coefficient
    D    <- ifelse(1:(nlayers-1) <= idzmix, Depi, Dhypo)
    #Dzoo <- ifelse(1:(nlayers-1) <= idzmix, D + DIncrZooEpi, D + DIncrZooHyp)
    #if(time > 315) {
    #  D    <- c(0, D * SF(depths) * 2, 0) #no diffusion at upper and lower bounderies out of the system
    #}
    #else {
      D    <- c(0, D , 0) #no diffusion at upper and lower bounderies out of the system
    #}
    #Dzoo <- c(0, Dzoo, 0)
    ## transport loop over all state variables
    for(i in 1:nstates){
      ## increase diffusion for zooplankton by a constant factor due to own mobility
      #if (i == 6) {
      #  Dori <- D
      #  D    <- Dzoo
      #}
      ## positions of the i-th state variable for all layers nlayers
      ##   (ReacTran format)
      id <- (i + (i-1)*(nlayers-1)) : ((i + (i-1)*(nlayers-1)) + nlayers - 1)
      ## positions of the i-th state variable for all layers nlayers
      ##   (SALMOformat)
      #id <- i + ((1:(nlayers-1))*nstates)
      tx <- x[id] # the i-th state variable for all layers
      ## flux of oxygen from atmosphere into first layer
      ## !!!           .................  !!! make this flexible
      fluxup <- ifelse(i == iO2, K2 * (o2sat(forcings[itemp]) - tx[1]) / dz[1], 0)
      ## increase vsink for cyanos in autumn
      #if (time > 315) {
        #vmat[,3]     <- vmat[,3] * 5
        #vmatsedi[,3] <- vmatsedi[,3] * 5
        #print(time)
        #vmat     <- vmat * 2
        #vmatsedi <- vmatsedi * 2
      #}
      ## transport of the i-th state variable
            dxx[id] <- tran.1D(tx, C.up = 0, C.down = 0, D = D, v = vmat[,i],
                         flux.up = fluxup, A = Area, dx = dz)$dC
      #if (i == 6) D <- Dori
    }

    ## sedimentation into sediment (only calculated when external sedimentation is switched on / internal sedimentation is switched off)
    if(cc["SF"] == 0){

      depthsvec <- rep(depths, nstates)
      sedicheck <- ifelse(depthsvec > zres, TRUE, FALSE)
      sedi <- as.vector(vmatsedi[2:nrow(vmatsedi),]) / rep(dz, nstates) *  # vsink / dz
        rep(forcings[isf + (0:(nlayers-1)) * ni], nstates) *            # SF
        (1 - rep(forcings[iaver + (0:(nlayers-1)) * ni], nstates)) *      # 1 - aver
        x * sedicheck                                     # tief > zres
      dxx <- dxx - sedi
    }
    ## return derivatives
    dxx
  })
}
