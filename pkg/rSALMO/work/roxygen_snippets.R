

#' Call SALMO Shared Library
#' 
#' This function calls the shared library of SALMO (version with macrophyte
#' coupling).
#' 
#' 
#' @param cfunc string, name of the C function to be called
#' @param nOfVar vector with number of variables
#' @param cc vector of constants for the model
#' @param pp matrix of phytoplankton parameters
#' @param uu input vector (environmental conditions)
#' @param xx state vector
#' @param pm parameters for macrophytes
#' @param mx states of macrophytes (ReacTran format)
#' @return list, containing the derivatives as first element
NULL





#' Hypsographic Functions from Tabular Interpolation
#' 
#' The main function creates a function list object from a data table to allow
#' interpolation of depths, areas and volume fractions of layers.
#' 
#' 
#' @param table hypsographic table (bathymetry) given as data frame with
#' columns \code{level}, \code{area} and \code{volume}
#' @return a list of functions with the following elements (slots):
#' 
#' \itemize{ \item \code{level(volume)} interpolates level from volume
#' 
#' \item \code{volume(level)} interpolates volume from a given water level
#' 
#' \item \code{area(level)} interpolates area from a given level
#' 
#' \item \code{vh(level, zmixreal)} calculates hypolimnic volume from surface
#' water level (m above ground) and mixing depth (zmixreal, thermocline in m
#' below surface)
#' 
#' \item \code{ve(level, zmixreal)} calculates epilimnic volume from surface
#' water level and zmixreal
#' 
#' \item \code{zmix(level, zmixreal)} calculates average mixing depth (zmix)
#' from surface water level and zmixreal
#' 
#' \item \code{zhm(level, zmixreal)} calculates average depth of bottom layer
#' (hypolimnion and metalimnion)
#' 
#' \item \code{sediment_area(level)} calculates sediment area for a series of
#' water levels, where \code{level} shoudl be a vecor of water levels in
#' ascending order
#' 
#' \item \code{pelagic_ratio(level)} calculates the area ratio of the pelagic
#' part (water) to the total area of layers. It is equal to 1 - sediment area
#' fraction.
#' 
#' }
#' @seealso \code{\link{bautzen_hypso}}, an example of a hypsographic table \cr
#' \code{\link{areaFunction}}, functions to create a hypsographic table for
#' lakes with canonical shape.
#' @examples
#' 
#' data(bautzen_hypso)                          # load data table
#' 
#' hypso <- set_hypso_functions(bautzen_hypso)  # create function object
#' 
#' hypso$volume(166)                # volume at 166 m a.s.l.
#' 
#' hypso$sediment_area(c(160, 166)) # sediment area of hypo- and epilimnion
#' 
NULL


