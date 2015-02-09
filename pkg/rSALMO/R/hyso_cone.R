#' Hypsographic Table of a Conic Lake
#'
#' Hypsographic Table of a Lake with Conical Morphometry
#' 
#' @source Generated data, see example below.
#' 
#' @references
#' 
#' Sachse, R., Petzoldt, T., Blumstock, M., Moreira, S., Pätzig, M., 
#'   Rücker, J., Janse, J. H., Mooij, W. M., Hilt, S. (2014)
#'   Extending one-dimensional models for deep lakes to simulate the impact of 
#'   submerged macrophytes on water quality.
#'   Environmental Modelling & Software, 61, 410-423
#'
#' @name hypso_cone
#' @docType data
#' @keywords data
#'
#' @examples
#' 
#' area   <- function(Amax, zmax, z){Amax * (1 - (z)/zmax)^2}
#' volume <- function(Amax, zmax, z){Amax * zmax/3 * (1 - z/zmax)^3}
#' zmax   <- 30                            # maximum depth (m)
#' Amax   <- 12.07 * 1000 * 1000           # total Area (m^2)
#' depth  <- seq(30, 0, -.1)               # depth (m)
#' 
#' hypso_cone <- data.frame(depth = depth, 
#'   volume = volume(Amax, zmax, depth), 
#'   area=area(Amax, zmax, depth)
#' )
#'
#' ## inverted  depth order
#' hypso_cone$level <- zmax - hypso_cone$depth
#' 
#' 
NULL