#' Column Names of SALMO State and Input Variables
#' 
#' 
#' 
#' @param macrophytes logical if state variable names for macrophytes should be included
#' @param nlayers numeric, number of layers fro replicating the names
#' @param parms list containing constant model parameters
#' 
#' @return vector of names for states resp inputs.
#' contain optional outputs
#'
#' @rdname salmo_state_names
#' @export salmo_input_names
#'
salmo_input_names <- function() {
  c("time", "vol", "depth", "dz", "qin", "ased", "srf",
  "iin", "temp", "nin", "pin", "pomin", "zin", "oin", "aver",
  "ad", "au", "eddy", "x1in", "x2in", "x3in", "sf")
}

#' @rdname salmo_state_names
#' @export salmo_state_names
#'
salmo_state_names <- function(nlayers, macrophytes = FALSE) {
  salmo_states <- c("N", "P", "X1", "X2", "X3", "Z", "D", "O", "G1", "G2", "G3")
  mac_states <- c("DVeg", "PVeg", "NVeg", "fRootVeg")
  if (macrophytes) salmo_states <- c(salmo_states, mac_states)
  rep(salmo_states, each = nlayers)
}
