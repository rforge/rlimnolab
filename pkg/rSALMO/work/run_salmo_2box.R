#void SalmoKern(int* nOfVar, double* c, double* p, double* u, double* x, double* dxq, double* dxs)
library(deSolve)
dyn.load("SalmoKern.dll")
source("SalmoNCP.R")

call_SalmoDLL <- function(cfunc,         # aufzurufene C-Funktion
                          nOfVar,        # Anzahl der Variablen
                          cc,            # Vektor der Modellkonstanten
                          pp,            # Matrix der Phytoplanktonkenngroessen
                          uu,            # Vektor der Randbedingungen
                          xx) {          # Zustandsvektor

  dxq <- dxs <- numeric(length(xx))
  ret  <- .C(cfunc, as.integer(nOfVar), cc=as.double(cc),
             pp=as.double(pp), uu=as.double(uu), xx=as.double(xx), 
             dxq=as.double(dxq), dxs=as.double(dxs))
  list(#ret$uu, 
       ret$xx, ret$dxq, ret$dxs)#, ret$pp)
}


SALMO <- function(t, x, p) {
  #cat(t, "\n")
  x1 <- x[1:nstates]
  x2 <- x[(nstates+1):(2*nstates)]
  ret1 <- call_SalmoDLL("SalmoCore", nOfVar, cc, pp, uu, x1)
  ## adapt depth, volume, surface area, "surface" light ...
  ret2 <- call_SalmoDLL("SalmoCore", nOfVar, cc, pp, uu, x2)
  ## return source - sink terms
  dx1 <-c(ret1[[2]] - ret1[[3]])
  dx2 <-c(ret2[[2]] - ret2[[3]])
  list(c(dx1, dx2))
}

## =============================================================================

nOfVar["numberOfLayers"]     <- 1
nOfVar["numberOfStates"]     <- 8 # ignored??
nOfVar["numberOfParameters"] <- length(pp) / nOfVar["numberOfAlgae"] ## important!
nOfVar["timestep"] <- 10

nstates <- nOfVar["numberOfStates"]

uu <- read.table("utest.txt", header=TRUE)
uu <- uu[1,]

uu["upomin"] <- 0
uu["upin"]   <- 20
#uu["uzmix"] <-  10  # ignored ?
#uu["utief"]  <- 1  # for resuspension

cc["EPSMIN"] <- 0.3 # ignored !? maybe in the physics module
cc["APSFMIN"] <- 1  # P release from sediment
cc["APSFMAX"] <- 7  # P release from sediment

#       N  P   X1   X2  X3 Z    D   O
x0 <- c(N=5, P=0, X1=1e-6, X2=1e-6,  X3=1e-6, Z=1e-6, D=0, O=12)
x0 <- c(N=5, P=10, X1=1e-6, X2=1e-6,  X3=1e-6, Z=1e-6, D=0, O=12)

x0 <- c(x0, x0)

## call one time for testing
ret <- call_SalmoDLL("SalmoCore", nOfVar, cc, pp, uu, x0)


times <- seq(0, 100, 1)

system.time(
  out <- ode(x0, times, SALMO, parms = NULL, method="lsoda")
)

plot(out)



















dyn.unload("SalmoKern.dll")