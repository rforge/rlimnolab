#' Hypsographic Functions from Tabular Interpolation
#' 
#' The main function creates a function list object from a data table
#' to allow interpolation of depths, areas and volume fractions of layers. 
#' 
#' @param table hypsographic table (bathymetry) given as data frame with columns \code{level}, 
#'        \code{area} and \code{volume} 
#'        
#' @return  a list of functions with the following elements (slots):
#' 
#' \itemize{
#'   \item \code{level(volume)} interpolates level from volume
#'   
#'   \item \code{volume(level)} interpolates volume from a given water level
#'   
#'   \item \code{area(level)} interpolates area from a given level
#'   
#'   \item \code{vh(level, zmixreal)} calculates hypolimnic volume from surface
#'         water level (m above ground) and mixing depth (zmixreal, 
#'         thermocline in m below surface)
#'         
#'   \item \code{ve(level, zmixreal)} calculates epilimnic volume from surface water 
#'         level and zmixreal
#'         
#'   \item \code{zmix(level, zmixreal)} calculates average mixing depth (zmix) from surface 
#'         water level and zmixreal
#'         
#'   \item \code{zhm(level, zmixreal)} calculates average depth of bottom layer 
#'         (hypolimnion and metalimnion)
#'         
#'   \item \code{sediment_area(level)} calculates sediment area for a series of water levels,
#'         where \code{level} shoudl be a vecor of water levels in ascending order
#'         
#'   \item \code{pelagic_ratio(level)} calculates the area ratio of the pelagic part (water)
#'        to the total area of layers. It is equal to 1 - sediment area fraction.
#' }
#' 
#' @seealso \code{\link{bautzen_hypso}}, an example of a hypsographic table \cr
#'          \code{\link{areaFunction}}, functions to create a hypsographic table 
#'             for lakes with conical shape.
#' 
#' @examples
#' 
#' data(bautzen_hypso)                      # load data table
#' hypso <- hypso_functions(bautzen_hypso)  # create function object
#' 
#' hypso$volume(166)                        # volume at 166 m a.s.l.
#' hypso$sediment_area(c(160, 166))         # sediment area of hypo- and epilimnion


hypso_functions <- function(table) {
  level  <- with(table, approxfun(volume, level))
  area   <- with(table, approxfun(level, area))
  volume <- with(table, approxfun(level, volume))
  
  vh <- function(level, zmixreal) {
    volume(level - zmixreal)
  }
  
  ve <- function(level, zmixreal) {
    vh <- vh(level, zmixreal)
    v  <- volume(level)
    v - vh
  }
  
  zmix <- function(level, zmixreal) {
    v <- volume(level)
    vh <- vh(level, zmixreal)
    a  <- area(level)
    (v - vh) / a
  }
  
  zhm <- function(level, zmixreal, ah.min = 0.1) {
    vh  <- vh(level, zmixreal)
    ah  <- area(level)
    if (ah > ah.min)
      zhm <- vh /ah
    else
      zhm <- 0
    zhm
  }
  
  ## sediment contact area of a layer; especially useful if level > 1
  sediment_area <- function(level) {
    if (length(level) < 1) stop("level is empty")
    a <- area(level)
    c(a[1], diff(a))
  }
  
  pelagic_ratio <- function(level){
    total_area    <- area(level)
    sediment_area <- c(total_area[1], diff(total_area))
    1 - sediment_area / total_area
  }
  
  list(level = level, area = area, volume = volume,
       ve = ve, vh = vh, zmix = zmix, zhm = zhm,
       sediment_area = sediment_area,
       pelagic_ratio = pelagic_ratio)
}







