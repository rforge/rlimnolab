### ============================================================================
### Demonstration of a 1box simulation
### Here we use the data set of bautzen reservoir, but simulate only
### an upper layer with fixed depth, just for testing.
### See 2 box version for a more realistic case.
### ============================================================================

### Warning: transport, sedimentation and sediment not yet implemented

library(rSALMO)

## Data set from workgroup limnology of TU Dresden
data(data_bautzen_1997)


## Reformat old-style input data structure into new structure 
forc <- with(data_bautzen_1997, 
             data.frame(
               time   = t,          # simulation time (in days)
               vol    = v,          # volume (m^3)
               depth  = s - 154,    # actual depth of the lake (m)
               dz     = 5,          # zmix, or layer depth (m)
               qin    = qin,        # water inflom (m^3 d^-1)
               ased   = 0,          # sediment contact area of the layer (m^2 ??)
               srf    = srf,        # strong rain factor, an empirical index of turbidity 
               iin    = iin,        # photosynthetic active radiation (J cm^2 d^-1); approx 50% of global irradiation
               temp   = temp,       # water temperature (deg. C)
               nin    = nin,        # DIN concentration in inflow (mg L^-1)   
               pin    = pin,        # DIP concentratioin in inflow (mug L^-1)
               pomin  = pomin,      # particulate organic matter in inflow, wet weight (mg L^-1)
               zin    = zin,        # zooplankton in inflow (w.w. mg L^-1)
               oin    = o2sat(temp),  # oxygen concentration in inflow (mg L^-1)
               aver   = 1,          # ratio of sediment contact area to total area
               ad     = 0,          # downwards flux between layers (m^3 d^-1)  
               au     = 0,          # upwards flux between layers (m^3 d^-1)
               diff   = 0,          # eddy diffusion coefficient 
               x1in   = 0, # xin1,       # phytoplankton import of group 1 (w.w. mg L^-1)
               x2in   = 0, # xin2,       # phytoplankton import of group 2 (w.w. mg L^-1)       
               x3in   = 0           # phytoplankton import of group 3 (w.w. mg L^-1)
             )
)
## Matrices are faster than data frames
forc <- as.matrix(forc)

## Read default parameter set. The order of the vector must not be changed!
parms <- get_salmo_parms()


## A few parameters that are specific for Bautzen Reservoir
parms$cc[c("MOMIN",  "MOT", "KANSF", "NDSMAX",	"NDSSTART",	"NDSEND",	"KNDS",	"KNDST")] <-
   c(0.005,   0.002,	    0,	 0.095,	   0,	         365,	     0.00,	 1.03)

## Background light extinction is lake specific
parms$cc["EPSMIN"] <- 0.7


## Initial values
## Xi = Phytoplankton biomass (mg/L w.w.)
## Z = Zooplankton Biomass    (mg/L w.w.)
## D = allochthonous detritus (mg/L w.w.)
## O = Oxygen concetration (mg/L)
## Gi = Carbon:Chlorophyll ratio for Baumert's photosynthesis model
##      (currently not used)
x0 <- c(N=5, P=10, X1=.1, X2=.1,  X3=.1, Z=.1, D=20, O=14, G1=0, G2=0, G3=0)


## Call one time for testing
ret <- call_salmodll("SalmoCore", parms$nOfVar, parms$cc, parms$pp, forc, x0)

## Simulation time steps
times <- seq(0, 365, 1)

## Model simulation with "lsoda" from package deSolve
#system.time(
  out <- ode(x0, times, salmo_1box, parms = parms, method="lsoda", inputs=forc)
#)

plot(out, which=c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "iin"))


