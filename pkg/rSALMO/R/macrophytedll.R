macrophytedll <-
function(
  cfunc,         # C function to be called
  mnOfVar,       # number of variables in Macrophyte Module
  nOfVar,        # number of variables in SALMO
  p,             # vector of macrophyte parameters
  cc,            # vector of macrophyte constants
  salmouu,       # vector of SALMO boundary conditions
  time,          # time [d]
  salmox,        # vector of SALMO states (all layers)
  mx             # vector of macrophyte states (all layers)
){
  ## create empty data structure for derivatives
  salmodmx <- numeric(length(salmox))
  dmx      <- numeric(length(mx))
  ## call C function
  ret <- .C(
    cfunc,
    mnOfVar = as.double(mnOfVar),
    nOfVar  = as.double(nOfVar),
    p       = as.double(p),
    cc      = as.double(cc),
    salmouu = as.double(salmouu),
    time    = as.double(time),
    salmox  = as.double(salmox),
    salmodmx= as.double(salmodmx),
    mx      = as.double(mx),
    dmx     = as.double(dmx)
  )
  ## return
  ## sort salmodmx for use with ReacTran 
  list(arrangeStateWise(ret$salmodmx, 
    nOfVar["numberOfStates"], nOfVar["numberOfLayers"]), ret$dmx)
}
