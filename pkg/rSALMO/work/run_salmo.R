library(marelac)
library(rSALMO)

data(bautzen1997)

data(cc)  # unique SALMO parameters, vector
data(pp)  # phytoplankton parameters, matrix with 1 column per phytopl. species 


forc <- with(bautzen1997, 
  data.frame(
    time   = t,
    vol    = ve, # v
    depth  = s - 154,
    dz     = 5,  # zmix,
    qin    = qin,
    ased   = ve / 5, # ve / zmix OR v/zmix
    srf    = 0,     # <-----
    iin    = iin,
    temp   = temp,
    nin    = nin,
    pin    = pin,
    pomin  = pomin,
    zin    = zin,
    oin    = o2sat(temp),
    aver   = 0,      # open water part of bottom area (1 - a_sed/a_layer) 
    ad     = 0,
    au     = 0,
    diff   = 0,
    x1in   = xin1,
    x2in   = xin2,
    x3in   = 0
  )
)
forc <- as.matrix(forc)
## =============================================================================

# NOFVAR
nOfVar <- c(
  numberOfInputs      = 21,
  numberOfOutputs     = 12,
  numberOfStates      = 11,
  numberOfParameters  = NA,
  numberOfAlgae       = 3,
  numberOfLayers      = 1,
  numberOfTributaries = 1,
  numberOfOutlets     = 1,
  timestep = 1
) 

nOfVar["numberOfParameters"] <- length(pp) / nOfVar["numberOfAlgae"] ## important!


cc["EPSMIN"] <- 0.7

cc[c("MOMIN",	"MOT", "KANSF", "NDSMAX",	"NDSSTART",	"NDSEND",	"KNDS",	"KNDST")] <-
   c(0.005,   0.002,	    0,	 0.095,	   0,	         365,	     0.00,	 1.03)

#       N  P   X1   X2  X3 Z    D   O
x0 <- c(N=5, P=10, X1=.1, X2=.1,  X3=.1, Z=.1, D=20, O=14, G1=0, G2=0, G3=0)


## call one time for testing
#ret <- call_SalmoDLL("SalmoCore", nOfVar, cc, pp, forc, x0)


times <- seq(0, 365, 1)

#system.time(
  out <- ode(x0, times, SALMO.1box, parms = NULL, method="lsoda", inputs=forc)
#)

plot(out, which=c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "iin"))

forc2 <- forc
forc2[, "pin"] <- 0.1 * forc[, "pin"]

out2 <- ode(x0, times, SALMO.1box, parms = NULL, method="lsoda", inputs=forc2)

plot(out, out2, which=c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "iin"))


## Todo: 
#  - aver, SF
#  - ICE


