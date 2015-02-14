library(rSALMO)
library(plot3D)

parms <- get_salmo_parms(nlayers=60, macrophytes=TRUE)

parms$pp_ma["kMortVegSum", 1]  <- 0.01 # [d^-1] Background mortality, little bit increased
parms$pp_ma["MaxWaveMort", 1]  <- 0    # [d^-1] MortWaveMax = 0  -> no wave mortality

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

dz    <- diff(depth)[1] # all dz must be equal; unequal depths not yet implemented
level <- maxdepth - depth + dz
vol   <- hyps$vol(level)

depthmatrix <- matrix(rep(depth, each=ntime), nrow=ntime)

# entrainment of tributaries; specific depths, Q dependend on residence time
qinmatrix   <- matrix(0, nrow=ntime, ncol=ndepth)
qinmatrix[,1:11] <- rep(vol[21:31]/300, each=ntime) 


## global irradiation data
data(irad)
## aggregate hourly data to daily sums (J/cm^2/d)
daily <- aggregate_daily(irad$time, irad$irad2, basedate= "1970-01-01 00:00.00 UTC")


## assumption: Jan + Feb with ice
iin <- daily$y * (1 - 0.9* (daily$day < as.POSIXct("2005-03-01")))


## sedimentation speed for each state in different depths
## ?? rename to advectionmatrix?

vmatsedi <- sedimentation_matrix(parms, nstates, nphy, depth)
vmatsedi[61, ] <- vmatsedi[60, ] ## fixme

## sediment area
asedmatrix  <- matrix(rep(hyps$sediment_area(level), each=ntime), nrow=ntime)

## pelagic ratio; todo: rename aver to something else; pelratio?
## !! note pmax 0.1
avermatrix  <- matrix(rep(pmax(0.1, hyps$pelagic_ratio(level)), each=ntime), nrow=ntime)
avermatrix[,60] <-0 # no water at bottom -- check this !!!


parms2 <- c(parms,
  list(
    K2 = 1.0,             # oxygen reaeration from atmosphere
    depths = depth,  # !!! ToDo: make this obsolete
    nstates = nstates,
    nlayers = nlayers,
    ni = ni,              # number of inputs
    zresmax = 20          # maximum resuspension depth
  )
)


inputs <- list(
  time  = time,
  vol   = vol,
  depth = depthmatrix, # vector "depth" should be enough here
  dz    = dz, #dzmatrix,    # scalar 0.5  
  qin   = qinmatrix,
  ased  = asedmatrix,
  srf   = 1,
  iin   = iin * 0.5,       # conversion to PAR
  temp  = turbulence$temp, # tempmatrix,
  nin   = 10,
  pin   = 50,
  pomin = 50,
  zin   = 1,
  oin   = 10,
  aver  = avermatrix,
  ad    = 0,
  au    = 0,
  eddy  = turbulence$eddy, #eddymatrix,
  x1in  = 5,
  x2in  = 5,
  x3in  = 5,
  sf    = 1, # no sediment focussing; note also vmatsedi
  vmatsedi = vmatsedi
)


## forcing **function**
forcing_functions <- function(inputs) {
  cnt <- 0

  ## otional:do some precalculations ...
  colnames <- salmo_input_names()
  
  forcings = function(time) {
    ## print simulation time to screen
    cnt <<- cnt + 1
    if (cnt > 100) {cnt <<- 0; cat("time := ", time, "\n")}

    ## ..... do interpolation here .....
    forc <- makeInputVector(inputs, time)
    
    attr(forc, "colnames") <- colnames
    
    ## check data and fix inconsistencies
    ## - temperature must be > 0
    forc[seq(9, 1320, 22)] <- pmax(forc[seq(9, 1320, 22)], 0.1)
    forc
  }
}
signal <- forcing_functions(inputs)

nspec <- parms2$nOfVar["numberOfStates"] + parms2$nOfVar_ma["numberOfStates"]  


cat(file="logfile.log")
syslog <- FALSE


## without macrophytes
y0 <- c(N=10, P=8, X1=0.1, X2=0.1,  X3=0.1, Z=0.1, D=1, O=13, G1=1, G2=1, G3=1)
y0 <- rep(y0, each = nlayers)
names(y0) <- NULL

state_names <- salmo_state_names(nlayers, macrophytes = FALSE)
nspec <- parms2$nstates

times <- 0:365 

ndx <- init_salmo_integers(parms)

# print(system.time(
# out   <- ode.1D(y = y0, times = times, func = salmo_1d, parms = parms2, 
#                 method = "bdf", nspec = nspec, atol = 1e-4, rtol=1e-4, hini = 0.1,
#                 inputs = inputs, forcingfun = signal, ndx=ndx
#          )
# ))

## -----------------------------------------------------------------------------

## start values with macropohytes
y0 <- c(N=10, P=8, X1=0.1, X2=0.1,  X3=0.1, Z=0.1, D=1, O=13, G1=1, G2=1, G3=1)

## Initial values for the macrophytes sub-model
xm <- c(
  sDVeg        = 1*50,       # [gDW/m^2] 
  sPVeg        = 0.002*50,   # [gP/m^2] 
  sNVeg        = 0.02*50,    # [gN/m^2] 
  afRootVeg    = 0.6         # [-]
)

y0 <- c(y0, xm)

nstates <- length(y0)

y0 <- rep(y0, each = nlayers)
names(y0) <- NULL

state_names <- salmo_state_names(nlayers, macrophytes = TRUE)

nspec <- parms2$nOfVar["numberOfStates"] + parms2$nOfVar_ma["numberOfStates"] 

print(system.time(
  out   <- ode.1D(y = y0, times = times, func = salmo_mac_1d, parms = parms2, 
                  method = "bdf", nspec = nspec, atol = 1e-4, rtol=1e-4, hini = 0.1,
                  inputs = inputs, forcingfun = signal, ndx=ndx
  )
))




#save.image("tmp_status.Rdata")
#write.table(out, "out.txt", col.names=NA, row.names=TRUE)

xmids <- depth - 0.5/2
names <- c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "G1", "G2", "G3")
names <- c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "G1", "G2", "G3", "DVeg", "PVeg", "NVeg", "fRootVeg")


### plot of state variables (outputs of SALMO)

image(out, grid = xmids, ylim=c(30, 0), main=names, add.contour=TRUE, which=c(1:8))

i <- 3:6  # state variables
image(out, grid = xmids, which=i, main=names[i], ylim=c(30,0), legend=TRUE)


i <- c(1,2,7,8)  # state variables
image(out, grid = xmids, which=i, main=names[i], ylim=c(30,0), legend=TRUE)


## P
image(out, grid = xmids, which=2, main=names[2], ylim=c(30,0), legend=TRUE)


## vegetation
i <- 12:15
image(out, grid = xmids, which=i, main=names[i], ylim=c(30,0), legend=TRUE)


