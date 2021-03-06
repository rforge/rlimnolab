### ============================================================================
### Demonstration of a 2box simulation
### Here we use the data set of Bautzen reservoir with
### flexible mixing depth (ZMIX)
### ============================================================================
library(rSALMO)

### ================ work in progress ======================
### Warning: transport, sedimentation and sediment not yet implemented
### ToDo:
###    - implement details of function SALMO.2box
### Open Tasks
###    - sediment area
###    - depth
###    - oxygen
### ========================================================


## Data set from workgroup limnology of TU Dresden
data(data_bautzen_1997)

## Hypsographic table of Bautzen reservoir
data(hypso_bautzen)

hyps <- hypso_functions(hypso_bautzen)

#bautzen1997$s

## sediment area of hypo- and epilimnion

## upper boundaries of hypo- and epilimnion
levels <- with(data_bautzen_1997, cbind(hypo = s - zmixreal, epi = s))

## calculate sediment area row-wise for pairs of hypo and epilimnion depths
## "1" means row wise, t() means that result needs to be transposed
## result: 1st colum: hypolimnic sediment area, 2nd: epilimnic sediment area
# ToDo: make this a function
sediment_area <- t(apply(levels, 1, hyps$sediment_area))

## same for pelagic ratio
pelagic_ratio <- t(apply(levels, 1, hyps$pelagic_ratio))

## check pelagic_ratio calculation
# areaE <- hyps$area(levels[,2])
# areaH <- hyps$area(levels[,1])
# 
# sedH <- areaH
# sedE <- areaE - sedH
# 
# ratioH <- 1 - sedH/areaH
# ratioE <- 1 - sedE/areaE


## Reformat old-style input data structure into new structure
forcE <- with(data_bautzen_1997,
             data.frame(
               time   = t,          # simulation time (in days)
               vol    = ve,          # volume (m^3)
               depth  = zmixreal,    # absolute depth of the layer (m), required for resuspension depth
               dz     = zmix,          # zmix, or layer thickness (m)
               qin    = qin,        # water inflow (m^3 d^-1)
               ased   = sediment_area[,2],          # sediment contact area of the layer (m^2); important
               srf    = srf,        # strong rain factor, an empirical index of turbidity
               iin    = iin,        # photosynthetic active radiation (J cm^2 d^-1); approx 50% of global irradiation
               temp   = temp,       # water temperature (deg. C)
               nin    = nin,        # DIN concentration in inflow (mg L^-1)
               pin    = pin,        # DIP concentratioin in inflow (mug L^-1)
               pomin  = pomin,      # particulate organic matter in inflow, wet weight (mg L^-1)
               zin    = zin,        # zooplankton in inflow (w.w. mg L^-1)
               oin    = oin,  # oxygen concentration in inflow (mg L^-1)
               aver   = pelagic_ratio[,2],      # 1 - (ratio of sediment contact area to total area); (redundant if SF=0)
               ad     = ah,          # downwards flux between layers (m^3 d^-1)
               au     = ae,          # upwards flux between layers (m^3 d^-1)
               diff   = 0,          # eddy diffusion coefficient, not used in 2box setting
               x1in   = xin1,       # phytoplankton import of group 1 (w.w. mg L^-1)
               x2in   = xin2,       # phytoplankton import of group 2 (w.w. mg L^-1)
               x3in   = 0,       # phytoplankton import of group 3 (w.w. mg L^-1)
               SF = 1
             )
)
forcH <- with(data_bautzen_1997,
              data.frame(
                time   = t,
                vol    = vh,
                depth  = s - 154,
                dz     = zhm,
                qin    = qhin,
                ased   = sediment_area[,1],
                srf    = srf,
                iin    = 0,            # internally calculated from epilimnion
                temp   = temph,
                nin    = nhin,
                pin    = phin,
                pomin  = pomhin,
                zin    = zhin,
                oin    = ohin,
                aver   = 0,            # is zero in hypolimnion
                ad     = ah,
                au     = ae,
                diff   = 0,
                x1in   = xhin1,
                x2in   = xhin2,
                x3in   = 0,
                SF = 1
              )
)


## assumption: Jan + Feb with ice
forcE$iin <- forcE$iin * (1 - 0.9* (forcE$time < 60))


## Matrices are faster than data frames
forc <- as.matrix(cbind(forcE, forcH))


## Read default parameter set.
parms <- get_salmo_parms()

## A few parameters that are specific for Bautzen Reservoir
parms$cc[c("MOMIN",  "MOT", "KANSF", "NDSMAX",	"NDSSTART",	"NDSEND",	"KNDS",	"KNDST")] <-
   c(0.005,   0.002,	    0,	 0.095,	   0,	         365,	     0.00,	 1.03)

## TEST TEST TEST
#parms$cc[c("MOMIN",  "MOT", "KANSF", "NDSMAX",  "NDSSTART",	"NDSEND",	"KNDS",	"KNDST")] <-
#   c(0.005,   0.002,	    0,	 0,	   0,	         0,	     0.00,	 1.03)


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

## 2 boxes -----------------

x0 <- rep(x0, 2)

## Call one time for testing
ret <- call_salmodll("SalmoCore", parms$nOfVar, parms$cc, parms$pp, forc, x0)

## Simulation time steps
times <- seq(0, 365, 1)

ndx <- init_salmo_integers(parms)

## Model simulation with package deSolve: try "lsoda", "adams", "vode", "ode45"
## rtol adapted to circumvent numerical problem at time = 132
system.time(
 out <- ode(x0, times, salmo_2box, parms = parms, method="adams", inputs=forc, ndx=ndx, atol=1e-6, rtol=1e-2)
)

plot(out, which=1:8)

plot(out, which=12:19)


## Todo:
#  - aver, SF
#  - ICE


