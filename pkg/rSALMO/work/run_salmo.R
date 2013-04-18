library(marelac)
library(simlake)

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
    time   = t,
    vol    = v,
    depth  = s - 154,
    dz     = 5, #zmix,
    qin    = qin,
    ased   = 0,
    srf    = 0,     # <-----
    iin    = iin,
    temp   = temp,
    nin    = nin,
    pin    = pin,
    pomin  = pomin,
    zin    = zin,
    oin    = o2sat(temp),
    aver   = 0,
    ad     = 0,
    au     = 0,
    diff   = 0,
    x1in   = xin1,
    x2in   = xin2,
    x3in   = 0,
    sf     = 1      # <-----
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

cc[c("MOMIN",	"MOT", "KANSF", "NDSMAX",	"NDSSTART",	"NDSEND",	"KNDS",	"KNDST")] <-
   c(0.005,   0.002,	    0,	 0.095,	   0,	         365,	     0.00,	 1.03)

#       N  P   X1   X2  X3 Z    D   O
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


