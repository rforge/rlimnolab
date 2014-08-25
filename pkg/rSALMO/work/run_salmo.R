library(marelac)
library(rSALMO)


parms <- get_salmo_parms()

parms$cc["EPSMIN"] <- 0.7

parms$cc[c("MOMIN",  "MOT", "KANSF", "NDSMAX",	"NDSSTART",	"NDSEND",	"KNDS",	"KNDST")] <-
  c(0.005,   0.002,	    0,	 0.095,	   0,	         365,	     0.00,	 1.03)


data(data_bautzen_1997)
forc <- with(data_bautzen_1997, 
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



#       N  P   X1   X2  X3 Z    D   O
x0 <- c(N=5, P=10, X1=.1, X2=.1,  X3=.1, Z=.1, D=20, O=14, G1=0, G2=0, G3=0)


## call dll function one time for testing
ret <- call_salmodll("SalmoCore", parms$nOfVar, parms$cc, parms$pp, forc, x0)

## test model function
salmo_1box(1, x0, p=parms, inputs=forc2)

times <- seq(0, 365, 1)



out <- ode(x0, times, salmo_1box, parms = parms, method="lsoda", inputs=forc)

plot(out, which=c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "iin"))

forc2 <- forc
forc2[, "pin"] <- 0.1 * forc[, "pin"]


out2 <- ode(x0, times, salmo_1box, parms = parms, method="lsoda", inputs=forc2)

plot(out, out2, which=c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "iin"))


## Todo: 
#  - aver, SF
#  - ICE


