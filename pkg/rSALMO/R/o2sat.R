#' Oxygen Saturation in Fresh Water
#'
#' Simple formula for oxygen saturation (in mg/L) from temperature (in deg C).
#'
#' @param T vector of temperatures (in deg C)
#'
#' @value oxygen saturation in mg/L
#'
#' @details This is a simple function can also be used to test availability
#'          of the SALMO shared library
#'
#' @seealso \code{\link[marelac]{gas_O2sat}} in package marelac
#'
#' @examples
#'
#' o2sat(c(10, 20))

o2sat <-
function(T) {
  sat1 <- function(T) {
    res <- numeric(length(T))
    ret <- .C("osat", T = as.double(T), result = as.double(res))
    ret$result
  }
  sapply(T, sat1)
}
