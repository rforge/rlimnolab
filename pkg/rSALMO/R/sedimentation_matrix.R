#' Create Sedimentation and Migration Matrix for SALMO-1D
#'
#' Note: state variable indexes still hard coded
#'
#' ToDo: pass parameters regularly to this function
#'
#' @param parms     SALMO parameter list
#' @param nstates   total number of state variables  (deprecated !!!)
#' @param nphy      number of phytoplankton groups   (deprecated !!!)
#' @param depth     vector of individual depths (m)  (deprecated !!!)
#' @param focussing optional two-valued vector for a linear sediment focussing 
#'   heuristics (e.g. \code{c(1, 2)} - experimental!)
#' @return matrix with sedimentation / migration velocities
#'
#' @examples
#' 
#' nstates <- 8 # nOfVar$numberOfStates
#' nphy    <- 3 # nOfVar$numberOfAlgae
#' parms   <- get_salmo_parms(nlayers=140, macrophytes=TRUE)
#' depth  <- seq(0, 70, 0.5) # unique(inputs[,3])
#' vmat   <- sedimentation_matrix(parms, nstates, nphy, depth)
#' 
#' @export sedimentation_matrix
#'

sedimentation_matrix <- function(parms, nstates=NULL, nphy=NULL, depth=NULL, focussing = c(1, 1)) {
  if (is.null(nstates)) nstates <- parms$nstates
  if (is.null(nphy))    nphy    <- parms$nphy
  if (is.null(depth))   depths  <- parms$depth
  
  nlayers <- length(depth)
  
  v <- numeric(nstates)
  v[3:(2+nphy)] <- parms$pp["VS",] # phytoplankton sedimentation velocities
  v[3+nphy+1]   <- parms$cc["VD"]  # detritus sedimentation velocity
  
  ## Sediment focussing heuristics: increasing sedimentation velocity with depth
  SF            <- approxfun(range(depth), focussing, rule = 2)
  vmat          <- matrix(rep(v, nlayers), byrow = TRUE, nrow = nlayers)
  vmat          <- vmat * SF(depth)
  vmat[,6]      <- 0 # -0.15 # - parms$cc["VMIG"] # measurement unit?
  ### add row with zeros at bottom
  #vmat          <- rbind(vmat, matrix(rep(0, nstates), byrow = TRUE, nrow = 1))
  ## or: additional element same as the one before
  vmat          <- rbind(vmat, vmat[nrow(vmat), ])

  ## no flux down of phytoplankton in last layer,
  ## because sedimentation to sediment is handled by SalmoCore internally
  #vmat[nrow(vmat),] <- rep(0,nstates)

  vmat
}


