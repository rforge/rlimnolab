library(rSALMO)

parms <- get_salmo_parms(nlayers=61, macrophytes=TRUE)

parms$pp_ma["kMortVegSum", 1]  <- 0.01 # [d^-1] Background mortality, little bit increased
parms$pp_ma["MaxWaveMort", 1]  <- 0    # [d^-1] MortWaveMax = 0  -> no wave mortality
#parms$pp_ma["pWaveMort", 1] <- 3      # [-] exponent for steepness of sigmoidal curve of wave mortality
#parms$pp_ma["HWaveMort" ,1] <- 0.5    # [m] Half saturation depth for wave mortality



parms$pp["EPSX",] <-  c(0.03, 0.02, 0.02) * 5 # specific extinction coefficients (from AEMON database, converted from mg Chl m^-3 to g WW m^-3 by factor of 5)
parms$pp["VS",] <- c(0.1, 0.4, 0.4)        # [m/d] settling velocity for phyto-plankton (?? from Rinke et. al. 2010)
parms$pp["VS",] <- c(0.1, 0.2, 0.2)        # [m/d] settling velocity for phyto-plankton

parms$cc["SEZMAX"] <- 0.05 # 0.4
parms$cc["APSFMIN"] <- 2
parms$cc["APSFMAX"] <- 7
parms$cc["EPSMIN"] <- 0.1

# SALMO calculation of sedimentation to sediment: internal  = 1, external = 0
parms$cc["SF"] <- 0

parms$nOfVar["numberOfInputs"] <- 22

nstates  <- parms$nOfVar["numberOfStates"]
nlayers  <- parms$nOfVar["numberOfLayers"]
nphy     <- parms$nOfVar["numberOfAlgae"]
ni       <- parms$nOfVar["numberOfInputs"]

## hypsographic function
data(hypso_cone)
hyps <- hypso_functions(hypso_cone)

###-----------------------------------------------------------------------------
### inputs
###-----------------------------------------------------------------------------
data(turbulence)

time        <- turbulence$time
depth       <- turbulence$depth
ntime       <- length(time)
ndepth      <- length(depth)
maxdepth    <- max(depth)

eddymatrix  <- turbulence$eddy * 8.64  # m2/s --> dm2/d to conversion
tempmatrix  <- turbulence$temp

depthmatrix <- matrix(rep(depth, each=ntime), nrow=ntime)
#emptymatrix <- matrix(0, nrow=ntime, ncol=ndepth)

dzmatrix    <- matrix(0.5, nrow=ntime, ncol=ndepth)
#qinmatrix   <- emptymatrix

## sediment area
asedmatrix  <- matrix(rep(hyps$sediment_area(depth), each=ntime), nrow=ntime)

## pelagic ratio; todo: rename aver to something else; pelratio?
avermatrix  <- matrix(rep(hyps$pelagic_ratio(depth), each=ntime), nrow=ntime)

### hier weiter ...

data(irad)

## derive daily sums (J/cm^2/d)
data(irad)
irad$day   <- floor(as.numeric(irad$time)/60/60/24)
daily      <- aggregate(list(irad = irad$irad2[-1]), list(day=irad$day[-1]), sum)
daily$time <- as.POSIXct((daily$day) * 60*60*24, origin = "1970-01-01 00:00.00 UTC")

daily <- data.frame(
  time = as.POSIXct(format(daily$time, "%Y-%m-%d")),
  irad = daily$irad
)

# ice and snow???
#iin <- forcings$iin * (1 - 0.9 * (icesnow$ice > 5))


## sedimentation speed for each state in different depths
## ?? advectionmatrix
vmatsedi <- sedimentation_matrix(parms, nstates, nphy, depth)


## start values for SALMO
y0 <- c(N=10, P=8, X1=0.1, X2=0.1,  X3=0.1, Z=0.1, D=1, O=13, G1=1, G2=1, G3=1)

## Initial values for the macrophytes sub-model
xm <- c(
  sDVeg        = 1*50,       # [gDW/m^2] 
  sPVeg        = 0.002*50,   # [gP/m^2] 
  sNVeg        = 0.02*50,    # [gN/m^2] 
  afRootVeg    = 0.6         # [-]
)

y0 <- rep(c(y0, xm), each = nlayers)
names(y0) <- NULL


parms2 <- c(parms,
  list(
    #Depi = 0,#50,  # zero or negative values switch "eddy" on
    #Dhypo = 0,# 1,  #
    K2 = 1.0,   # reaeration from atmosphere
    #cc = cc, pp = pp, 
    #cc_ma = cc_ma, pp_ma = pp_ma,
    #nOfVar = nOfVar, nOfVar_ma = nOfVar_ma,
    depths = depth,
    nstates = nstates,
    nlayers = nlayers,
    ni = ni   # number of inputs
  )
)

inputs <- list(
  time = time,
  vol  = hyps$vol(30),
  depth = depthmatrix, # vector would be enough here
  dz    = dzmatrix,    # scalar 0.5  
  qin = 0,
  ased=asedmatrix,
  srf=1,
  iin = daily$irad,
  temp = tempmatrix,
  nin =0,
  pin=0,
  pomin=0,
  zin=0,
  oin=0,
  aver = avermatrix,
  ad=0,
  au=0,
  eddy = eddymatrix,
  x1in=0,
  x2in=0,
  x3in=0,
  sf   = 1 + 4 * depthmatrix/maxdepth, # 1 ... 5
  vmatsedi = vmatsedi
)

tmp <- NULL

## forcing **function**
forcing_functions <- function(inputs) {
  cnt <- 0
  
  ## do some precalculations ...
  
  colnames <- salmo_input_names()
  
  forcings = function(time) {
    ## print simulation time to screen
    cnt <<- cnt + 1
    if (cnt > 100) {
      cnt <<- 0
      cat("time := ", time, "\n")
    }
    ## ..... do interpolation here .....
    forc <- makeInputVector(inputs, time)
    
    tmp <<- forc
    
    #forc <- makeInputVector(inputs, 0)
    
    attr(forc, "colnames") <- colnames
    
    ## check data and fix inconsistencies
    ## temperature must not be <= 0
    forc[seq(9, 1342, 22)] <- pmax(forc[seq(9, 1342, 22)], 0.1)
    forc
  }
}
signal <- forcing_functions(inputs)


### ..........
### ..........
### ..........

nspec <- parms2$nOfVar["numberOfStates"] + parms2$nOfVar_ma["numberOfStates"]  

#rm(cc, pp, cc_ma, pp_ma, nOfVar, nOfVar_ma)

times <- 0:365 # 365


cat(file="logfile.log")
syslog <- FALSE


## without macrophytes
y0 <- c(N=10, P=8, X1=0.1, X2=0.1,  X3=0.1, Z=0.1, D=1, O=13, G1=1, G2=1, G3=1)
y0 <- rep(y0, each = nlayers)
names(y0) <- NULL

state_names <- salmo_state_names(nlayers, macrophytes = FALSE)
nspec <- parms2$nstates

#tst <- salmo_1d(1, y0, parms2, inputs, forcingfun = signal)

out   <- ode.1D(y = y0, times = times, func = salmo_1d, parms = parms2, 
                method = "bdf", nspec = nspec, atol = 1e-4, rtol=1e-4, hini = 0.1,
                inputs = inputs, forcingfun = signal
         )



save.image("tmp_status.Rdata")
write.table(out, "out.txt", col.names=NA, row.names=TRUE)

xmids <- 30-depth - 0.5/2
names <- c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "G1", "G2", "G3")
#names <- c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "DVeg", "PVeg", "NVeg", "fRootVeg")


### plot of state variables (outputs of SALMO)

image(out, grid = xmids, ylim=c(30, 0), main=names, add.contour=TRUE, which=c(1:8))

i <- 3:5  # state variables
j <- 1    # which year
n <- 365  # length of year
image(out, grid = xmids, which=i, main=names[i], ylim=c(30,0), legend=TRUE, subset = ((j-1) * n + 1) : (j*n))

## reorder layer wise if you want ...
#ii <- as.vector(outer((0:139)*11, 1:11, "+")) + 1
#oo <- out[,ii]



### plot of SALMO inputs
#library(OceanView)
library(plot3D)

image2D(inputs$temp)
image2D(inputs$depth)
image2D(inputs$ased)
image2D(inputs$aver)


