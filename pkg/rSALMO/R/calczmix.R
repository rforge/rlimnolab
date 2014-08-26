#' Calculate Mixing Depth
#'
#' Calculate mixing depth from temperature profile using a simple temperature threshold.
#'
#' @param temp vector of temperatures (from top to bottom (?))
#' @param depths vector of depths
#' @return list with id and depth of zmix
#'
#' @export calczmix

calczmix <- function(temp, depths) {
  nlayers    <- length(temp)
  t_tresh    <- 1 * (depths[2] - depths[1])
  idzmix     <- which(abs(diff(temp)) >= t_tresh)[1]
  if(is.na(idzmix)) idzmix <- nlayers
  zres       <- depths[idzmix]
  ## return values
  list(idzmix = idzmix, zres = zres)
}
