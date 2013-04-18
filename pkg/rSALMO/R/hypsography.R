#' Hypsographic Functions for Conical Lake Geometry
#' 
#' The functions calculate area resp. volume of a lake with an ellipsoidic depth profile.
#' 
#' @param level vector of water level above ground
#' @param zmax  maximum lake depth
#' @param amax  surface area of the lake surface for \code{level = zmax}
#' @param vol   vector of given volume(s)
#' 
#' @examples
#' 
#' level <- seq(0, 50, 1)
#' amax  <- 9e6
#' zmax  <- 50
#' area <- areaFunction(level, zmax, amax)
#' vol  <- volumeFunction(level, zmax, amax)
#' plot(vol, cumsum(area), type="l")
#' (vol - cumsum(area)) / vol

#' @rdname areaFunction
areaFunction   <- function(level, zmax, amax) amax * (1 - (zmax - level)/zmax)^2

#' @rdname areaFunction
volumeFunction <- function(level, zmax, amax) amax * level^3 / (3 * zmax^2)

#' @rdname areaFunction
vol2level <- function(vol, zmax, amax) (3 * vol * zmax^2 / amax)^(1/3)


