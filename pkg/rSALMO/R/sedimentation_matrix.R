#' Create Sedimentation and Migration Matrix for SALMO-1D
#'
#' Note: state variable indexes still hard coded
#'
#' ToDo: pass parameters regularly to this function
#'
#' @param nstates   total number of state variables
#' @param nphy      number of phytoplankton groups
#' @param depths    ????
#' @param focussing two-valued vector for a linear sediment focussing heuristics
#' @return matrix with sedimentation / migration velocities
#'
#' @example
#' 
#' nstates <- 8 # nOfVar$numberOfStates
#' nphy    <- 3 # nOfVar$numberOfAlgae
#' depths <- seq(0, 70, 0.5) # unique(inputs[,3])
#' vmat <- sedimentation_matrix(nstates, nphy, depths)
#' 


sedimentation_matrix <- function(nstates, nphy, depths, focussing = c(1, 5)) {
  nlayers <- length(depths)
  v <- numeric(nstates)
  v[3:(2+nphy)] <- pp["VS",] # phytoplankton sinking velocities
  v[3+nphy+1]   <- cc["VD"]  # detritus sinking velocity
  # Sediment Focussing: increasing sedimentation velocities with depth
  SF            <- approxfun(range(depths), focussing, rule = 2)
  vmat          <- matrix(rep(v, nlayers), byrow = TRUE, nrow = nlayers)
  vmat          <- vmat * SF(depths)
  vmat[,6]      <- 0 # -0.15 # -cc["VMIG"] # measurement unit?
  ### add row with zeros on top
  vmat          <- rbind(matrix(rep(0, nstates), byrow = TRUE, nrow = 1), vmat)

  ## no flux.down of phytoplankton in last layer,
  ## because sedimentation to sediment is handled by SalmoCore internally
  #vmat[nrow(vmat),] <- rep(0,nstates)

  vmat
}


