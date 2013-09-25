### ============================================================================
### Demonstration of a 1box simulation
### Here we use the data set of bautzen reservoir, but simulate only
### an upper layer with fixed depth, just for testing.
### See 2 box version for a more realistic case.
### ============================================================================

### Warning: transport, sedimentation and sediment not yet implemented

library(rSALMO)

## Data set from workgroup limnology of TU Dresden
data(bautzen1997)

## Reformat old-style input data structure into new structure 
forc <- with(bautzen1997, 
             data.frame(
               time   = t,          # simulation time (in days)
               vol    = v,          # volume (m^3)
               depth  = s - 154,    # actual depth of the lake (m)
               dz     = 5,          # zmix, or layer depth (m)
               qin    = qin,        # water inflom (m^3 d^-1)
               ased   = v/5,        # sediment contact area of the layer (m^2 ??)
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



## Read default parameter set. The order of the vector must not be changed!
data(pp)  # Phytoplankton parameters of SALMO, matrix with 1 column per phytopl. species 
data(cc)  # other SALMO parameters, vector
data(pp_ma) # macrophyte parameters
data(cc_ma) # other macrophyte parameters

## Matrices are faster than data frames
forc <- as.matrix(forc)
vsedi <- matrix(c(0,0,pp["VS",],0,0,0,0,0,0),nrow=1)
inputs <- list(forcings=forc, vsedi=vsedi)


## specify macrophyte heights
## by setting fDepth1Veg=0.5 and fDepth2Veg=1 macrophytes grow in the lower half of the box. By 
## setting fDepth1Veg=0 and fDepth2Veg=0.5 macrophytes grow only in the upper half of the box.
## For this example I assume that macrophytes grow in the lower 4.5m of the box.
## I need to let them grow so big, otherwise due to light limitation no macrophytes
## will grow in eutrophic Bautzen reservoir.
pp_ma["cMaxHeightVeg", 1] <- 4.5 # maximum height of macrophytes [m]
pp_ma["fDepth1Veg",1] <- 5/50    # upper fraction of the box up to which submerged macrophytes grow [-]
pp_ma["fDepth2Veg",1] <- 50/50  # lower fraction of the box down to which submerged macrophytes grow [-]


## A few parameters that are specific for Bautzen Reservoir
cc[c("MOMIN",  "MOT", "KANSF", "NDSMAX",	"NDSSTART",	"NDSEND",	"KNDS",	"KNDST")] <-
   c(0.005,   0.002,	    0,	 0.095,	   0,	         365,	     0.00,	 1.03)

## Background light extinction is lake specific
cc["EPSMIN"] <- 0.7



## Set of technical parameters, a.o. for dimensioning dynamic variables
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

## Total number of parameters passed to the model (very important)
nOfVar["numberOfParameters"] <- length(pp) / nOfVar["numberOfAlgae"]


## set of technical parameters for macrophyte model
data(nOfVar_ma)

## total number of layers 
nOfVar_ma["numberOfLayers"] <- nOfVar["numberOfLayers"]


## Put all 3 kinds of model parameters in a list
parms <- list(pp=pp, cc=cc, nOfVar=nOfVar, pp_ma=pp_ma, cc_ma=cc_ma)
## -----------------------------------------------------------------------------


## Initial values
## Xi = Phytoplankton biomass (mg/L w.w.)
## Z = Zooplankton Biomass    (mg/L w.w.)
## D = allochthonous detritus (mg/L w.w.)
## O = Oxygen concetration (mg/L)
## Gi = Carbon:Chlorophyll ratio for Baumert's photosynthesis model
##      (currently not used)
x0 <- c(
  N=5,   #[mg/l]
  P=10,  #[micro gram / l]
  X1=.1, #[mg WW/l]
  X2=.1, #[mg WW/l]
  X3=.1, #[mg WW/l]
  Z=.1,  #[mg WW/l]
  D=20,  #[mg WW/l]
  O=14,  #[mg/l]
  G1=0, 
  G2=0, 
  G3=0, 
  sDVeg        = 1*50,       # [gDW/m^2] 
  sPVeg        = 0.002*50,   # [gP/m^2] 
  sNVeg        = 0.02*50,    # [gN/m^2] 
  afRootVeg    = 0.6         # [-]
)


## Call one time for testing
ret <- call_salmodll("SalmoCore", nOfVar, cc, pp, forc, x0)

## Simulation time steps
times <- seq(0, 365, 1)

## Model simulation with "lsoda" from package deSolve
#system.time(
  out <- ode(x0, times, salmo_mac_1box, parms = parms, method="lsoda", inputs=inputs)
#)

plot(out, which=c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "sDVeg", "sPVeg", "sNVeg", "afRootVeg"))


