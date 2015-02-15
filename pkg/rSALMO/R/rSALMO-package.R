#' Ecological Simulation of Lakes
#'
#' \tabular{ll}{
#' Package:  \tab rSALMO\cr
#' Type:     \tab Package\cr
#' Version:  \tab 0.1\cr
#' Date:     \tab 2013-04-27\cr
#' License:  \tab  GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' An R implementation of the SALMO lake model\cr 
#' (SALMO: Simulation by means of An analytical Lake MOdel)
#' @name rSALMO-package
#' @aliases rSALMO rSALMO-package
#' @docType package
#' @author Rene Sachse and Thomas Petzoldt (package)\cr
#' 
#' Susanne Rolinski (C code of SALMO)\cr
#' 
#' Juergen Benndorf + (system of equations of SALMO)\cr
#' 
#' Maintainer: Thomas Petzoldt <thomas.petzoldt@@tu-dresden.de>\cr
#' @seealso \code{\link[deSolve:deSolve-package]{deSolve}}
#' @references see \url{http://hhbio.wasser.tu-dresden.de/projects/salmo/}
#' @keywords package
#' @examples
#' 
#' demo("demo_salmo_1box")
#' demo("demo_salmo_2box")
#' 
#' @useDynLib rSALMO
#' 
#' @import deSolve ReacTran plot3D RColorBrewer
NULL
