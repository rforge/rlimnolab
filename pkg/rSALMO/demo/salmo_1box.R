library(marelac)
library(rSALMO)

call_SalmoDLL <- function(cfunc,         # C function to be called
                          nOfVar,        # numbers of variables
                          cc,            # vector of parameters
                          pp,            # matrix of phytoplankton parameters
                          uu,            # vector of forcing functions
                          xx) {          # vector of initial states

  dxq <- dxs <- numeric(length(xx))
  ret  <- .C(cfunc, as.integer(nOfVar), cc=as.double(cc),
             pp=as.double(pp), uu=as.double(uu), xx=as.double(xx), 
             dxq=as.double(dxq), dxs=as.double(dxs))
  list(ret$xx, ret$dxq, ret$dxs)
}


SALMO.1box <- function(time, x, p, inputs) {
  #cat(time, "\n")

  ## interpolate data (this is the slowest part of the simulation)
  uu <- approxTime1(inputs, time)
  
  ## call the C model core 
  ret <- call_SalmoDLL("SalmoCore", nOfVar, cc, pp, uu, x)
  
  ## return source - sink terms
  dx <-c(ret[[2]] - ret[[3]])
  
  ## fix oxygen balance at surface to saturated value
  dx[8] <-  o2sat(uu["temp"]) - x[8]
  #list(dx, dO = unname(dx[8]))
  list(dx, iin = unname(uu["iin"]))
}

data(bautzen1997)

data(cc)  # unique SALMO parameters, vector
data(pp)  # phytoplankton parameters, matrix with 1 column per phytopl. species 


forc <- with(bautzen1997, 
  data.frame(
    time   = t,          # simulation time (in days)
    vol    = v,          # volume (m^3)
    depth  = s - 154,    # actual depth of the lake (m)
    dz     = 5,          # zmix, or layer depth (m)
    qin    = qin,        # water inflom (m^3 d^-1)
    ased   = 0,          # sediment contact area of the layer (m^2 ??)
    srf    = 0,          # strong rain factor, an empirical index of turbidity 
    iin    = iin,        # photosynthetic active radiation (J cm^2 d^-1); approx 50% of global irradiation
    temp   = temp,       # water temperature (deg. C)
    nin    = nin,        # DIN concentration in inflow (mg L^-1)   
    pin    = pin,        # DIP concentratioin in inflow (mg? L^-1)
    pomin  = pomin,      # particulate organic matter in inflow, wet weight (mg L^-1)
    zin    = zin,        # zooplankton in inflow (w.w. mg L^-1)
    oin    = o2sat(temp),  # oxygen concentration in inflow (mg L^-1)
    aver   = 0,          # ratio of sediment contact area to total area
    ad     = 0,          # downwards flux between layers (m^3 d^-1)  
    au     = 0,          # upwards flux between layers (m^3 d^-1)
    diff   = 0,          # eddy diffusion coefficient 
    x1in   = xin1,       # phytoplankton import of group 1 (w.w. mg L^-1)
    x2in   = xin2,       # phytoplankton import of group 2 (w.w. mg L^-1)       
    x3in   = 0           # phytoplankton import of group 3 (w.w. mg L^-1)
  )
)
forc <- as.matrix(forc)
## =============================================================================

# NOFVAR
nOfVar <- c(
  numberOfInputs      = 21,
  numberOfOutputs     = 14,
  numberOfStates      = 11,
  numberOfParameters  = NA,
  numberOfAlgae       = 3,
  numberOfLayers      = 1,
  numberOfTributaries = 4,
  numberOfOutlets     = 3,
  timestep = 1
) 

nOfVar["numberOfParameters"] <- length(pp) / nOfVar["numberOfAlgae"] ## important!


cc["EPSMIN"] <- 0.7

## a few parameters that are specific for Bautzen Reservoir
cc[c("MOMIN",	"MOT", "KANSF", "NDSMAX",	"NDSSTART",	"NDSEND",	"KNDS",	"KNDST")] <-
   c(0.005,   0.002,	    0,	 0.095,	   0,	         365,	     0.00,	 1.03)

## Initial values
## X = Phytoplankton biomass, Z = Zooplankton Biomass
## 
x0 <- c(N=5, P=10, X1=.1, X2=.1,  X3=.1, Z=.1, D=20, O=14, G1=0, G2=0, G3=0)


## call one time for testing
ret <- call_SalmoDLL("SalmoCore", nOfVar, cc, pp, forc, x0)


times <- seq(0, 365, 1)

#system.time(
  out <- ode(x0, times, SALMO.1box, parms = NULL, method="lsoda", inputs=forc)
#)

plot(out, which=c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "iin"))


## Todo: 
#  - aver, SF
#  - ICE


