

library(deSolve)                       # numerical integrator
library(ReacTran)                      # 1D transport
library(rSALMO)

data(cc)  # unique SALMO parameters, vector
data(pp)  # phytoplankton parameters, matrix with 1 column per phytopl. species 

data(cc_ma)
data(pp_ma)

pp_ma[32,1]  <- 0.01 # [d^-1] Background mortality, little bit increased
#pp_ma[40,1] <- 3    # [-] exponent for steepness of sigmoidal curve of wave mortality
#pp_ma[41,1] <- 0.5  # [m] Half saturation depth for wave mortality
pp_ma[42,1]  <- 0    # [d^-1] MortWaveMax = 0  -> no wave mortality


## simulation settings
data(nOfVar)
data(nOfVar_ma)

pp[2,] <-  c(0.03, 0.02, 0.02) * 5 # specific extinction coefficients (from AEMON database, converted from mg Chl m^-3 to g WW m^-3 by factor of 5)
pp[20,] <- c(0.1, 0.4, 0.4)        # [m/d] settling velocity for phyto-plankton (from Rinke et. al. 2010)

# SALMO internal sedimentation to sediment = 1, external sedimentation to sediment = 0
# don't use internal sedimentation anymore, 
# since increase of vsink with depth is not considered internally at the moment
# to make this possible again, one would need to change vsink parameters for X1..3
# for every call of a SALMO Layer in MReaktion()
cc["SF"] <- 0 


nOfVar["numberOfLayers"]       <- 140
nOfVar["numberOfStates"]       <- 11 #11# thpe: increased from 8 to 11 because of G
nOfVar["numberOfParameters"]   <- length(pp) / nOfVar["numberOfAlgae"]
nstates                        <- nOfVar["numberOfStates"]
nlayers                        <- nOfVar["numberOfLayers"]
nphy                           <- nOfVar["numberOfAlgae"]
ni                             <- nOfVar["numberOfInputs"]
nOfVar_ma["numberOfParameters"]<- 42

###-----------------------------------------------------------------------------
### inputs
###-----------------------------------------------------------------------------

### mixed Lake Stechlin 2009 / Lake Scharmuetzel test set
data(inputs)



## phytoplankton import
inputs[,20] <- 20  # X2
inputs[,21] <- 40  # X3


## ice or other light reflection in winter
inputs[,8] <- ifelse (inputs[,1] < 50 | inputs[,1] > 310, inputs[,8]*0.05, inputs[,8]) 

## check names: what is ad, au, **aver** ?
colnames(inputs) <- c("time", "vol", "depth", "dz", "qin", "ased", "srf", 
  "iin", "temp", "nin", "pin", "pomin", "zin", "oin", "aver", 
  "ad", "au", "eddy", "x1in", "x2in", "x3in", "sf")


## rearrange layers to wide format
alltimes <- unique(inputs[,1])
forcings <- matrix(0, nrow = length(alltimes), ncol = ncol(inputs) * nlayers)
for(i in 1:length(alltimes)){
  forcings[i, ] <- layers2vector(inputs[inputs[ ,1] == alltimes[i], ])
}


## attach the names of forcings as an attribut
attr(forcings, "colnames") <- colnames(inputs)

## sedimentation speed for each state in different depths
v <- numeric(nstates)
v[3:(2+nphy)]     <- pp[20,]  # phytoplankton sinking velocities
v[3+nphy+1]       <- cc["VD"] # detritus sinking velocity
depths            <- unique(inputs[,3])
SF                <- approxfun(range(inputs[,3]),c(5,10)) #increasing sedimentation velocities with depth
vmat              <- matrix(rep(v, nlayers), byrow=TRUE, nrow=nlayers)    
vmat              <- vmat*SF(depths)
vmat[,6]          <- 0#-0.15 #VMIG
vmat              <- rbind(matrix(rep(0,nstates), byrow=TRUE, nrow=1), vmat)
vmatsedi          <- vmat
vmatsedi[,6]      <- 0
vmat[nrow(vmat),] <- rep(0,nstates) # no flux.down of phytoplankton in last layer, because sedimentation to sediment is dealt by the SalmoCore itself

## save original vsink distribution 
# thpe: redundant
#vmatsediori       <- vmatsedi
#vmatori           <- vmat


## care about ice-cover and increase light reflection
## maybe use simple empirical ice model for future projections


######
### inits
######
##      N  P   X1   X2  X3 Z    D   O
#x0 <- c(N=0.4, P=8, X1=1e-6, X2=1e-6,  X3=1e-6, Z=0.1, D=0, O=13)

## thpe
y0 <- c(N=0.4, P=8, X1=1e-6, X2=1e-6,  X3=1e-6, Z=0.1, D=0, O=13, G1=1, G2=1, G3=1)


y0 <- rep(y0, nlayers)
y0 <- arrangeStateWise(y0, nOfVar["numberOfStates"], nOfVar["numberOfLayers"])

xm <- c(    
  sDVeg        = 1*50,       # [gDW/m^2] 
  sPVeg        = 0.002*50,   # [gP/m^2] 
  sNVeg        = 0.02*50,    # [gN/m^2] 
  afRootVeg    = 0.6         # [-]
)
xm <- rep(xm, each=nlayers)

y0 <- c(y0, xm)
names(y0) <- NULL

######
### simulation
######
## observed time steps
#times <- seq(0, max(thetime), 1)
times <- seq(100, 150, 1)


parms <- list(
  Depi = 50,
  Dhypo = 1,
  K2 = 1,
  cc = cc, pp = pp, 
  cc_ma = cc_ma, pp_ma = pp_ma,
  nOfVar = nOfVar, nOfVar_ma = nOfVar_ma,
  depths = depths,
  nstates = nstates,
  nlayers = nlayers,
  ni = ni
)

inputs <- list(
  vmat = vmat,
  vmatsedi = vmatsedi,
  forcings = forcings
)  


nspec <- nOfVar["numberOfStates"] + nOfVar_ma["numberOfStates"]  

rm(cc, pp, cc_ma, pp_ma, nOfVar, nOfVar_ma, nstates, nlayers)

# unew is a global variable used in "getuvec"
# !!! nlayers, nOfVar are global variables used in "transport"
# depths required below

tt <- system.time(
  out   <- ode.1D(y=y0, times=times, func=salmo_mac_1d, parms = parms, method="bdf", 
    nspec = nspec, inputs = inputs)
)
print(tt)

## plot results
#plot(out)

xmids <- depths - 0.5/2
#names <- c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "G1", "G2", "G3")
names <- c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "DVeg", "PVeg", "NVeg", "fRootVeg")
#pdf("../Stechlin/Stechlin-full-coupled-model.pdf", width=20, height=15)
image(out, grid = xmids, ylim=c(70, 0), main=names, add.contour=TRUE, which=c(1:8,12:15))
#dev.off()


#write.table(out, "out.csv", col.names=NA, row.names=TRUE, sep=";", dec=",")







