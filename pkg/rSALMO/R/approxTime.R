#' Linear Interpolation with Complete Matrices or Data Frames
#' 
#' 
#' Return a data frame, matrix or vector which linearly interpolates data from
#' a given matrix or data frame.
#' 
#' The functions can be used for linear interpolation with a complete matrix or
#' data frame. This can be used for example in the main function of an
#' \code{odeModel} to get input values at a specified time \code{xout}.
#' Versions \code{approxTime1} and \code{approxTimeEq} are less flexible 
#' (only one single value for \code{xout} resp. equidistant values of of time in \code{x}
#' and only linear interpolation) but have increased performance.  All interpolation
#' functions are faster if \code{x} is a matrix instead of a data frame.
#' 
#' @aliases approxTime approxTime1 approxTimeEq
#' @param x a matrix or data frame with numerical values giving coordinates of
#' points to be interpolated. The first column needs to be in ascending order
#' and is interpreted as independent variable (e.g. time), the remaining
#' columns are used as dependent variables.
#' @param xout a vector (or single value for \code{approxTime1}) of independend
#' values specifying where interpolation has to be done.
#' @param rule an integer describing how interpolation is to take place outside
#' the interval [min(x), max(x)]. If \code{rule} is 1 then \code{NA}s are
#' returned for such points and if it is 2, the value at the closest data
#' extreme is used.
#' @param ... optional parameters passed to \code{approx}.
#' @return \code{approxTime} returns a matrix resp. data frame of the same
#' structure as \code{x} containing data which interpolate the given data with
#' respect to \code{xout}.  \code{approxTime1} is a performance optimized
#' special version with less options than the original \code{approx} function.
#' It returns an interpolated vector.
#' @seealso \code{\link[stats]{approxfun}}
#' @keywords arith
#' @examples
#' 
#' inputs <- data.frame(time = 1:10, y1 = rnorm(10), y2 = rnorm(10, mean = 50))
#' input  <- approxTime(inputs, c(2.5, 3), rule = 2)
#' 
#' @export approxTime
approxTime <- function(x, xout, ...) {
  if (is.data.frame(x)) {x <- as.matrix(x); wasdf <- TRUE} else wasdf <- FALSE
  if (!is.matrix(x)) stop("x must be a matrix or data frame")
  m <- ncol(x)
  y <- matrix(0, nrow=length(xout), ncol=m)
  y[,1] <- xout
  for (i in 2:m) {
    y[,i] <- as.vector(approx(x[,1], x[,i], xout, ...)$y)
  }
  if (wasdf) y <- as.data.frame(y)
  names(y) <- dimnames(x)[[2]]
  y
}

#' @rdname approxTime
#' @export approxTime1
approxTime1 <- function (x, xout, rule = 1) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if ((!is.numeric(xout)) | (length(xout) != 1))
    stop("xout must be a scalar numeric value")
  if ((!is.numeric(rule)) | (length(rule) != 1))
    stop("rule must be a scalar numeric value")

  n <- nrow(x)
  if (xout >= x[n, 1]) {
    y <- c(xout, x[n, -1])
    if (rule == 1 & (xout > x[n + 1]))
      y[2:length(y)] <- NA
  }
  else if (xout <= x[1, 1]) {
    y <- c(xout, x[1, -1])
    if (rule == 1 & (xout < x[1]))
      y[2:length(y)] <- NA
  }
  else {
    i <- which.max(x[, 1] > xout)
    x1 <- x[i - 1, 1]
    x2 <- x[i, 1]
    y1 <- x[i - 1, ]
    y2 <- x[i, ]
    y <- y1 + (y2 - y1) * (xout - x1)/(x2 - x1)
  }
  names(y) <- dimnames(x)[[2]]
  y
}

#' @rdname approxTime
#' @export approxTimeEq
approxTimeEq <- function(x, xout, ...) {
  if (is.data.frame(x)) {x <- as.matrix(x); wasdf <- TRUE} else wasdf <- FALSE
  if (!is.matrix(x)) stop("x must be a matrix or data frame")
  n <- nrow(x)
  m <- ncol(x)
  nout <- length(xout)
  rxo <- range(xout, finite = TRUE)
  y <- matrix(NA, nrow = nout, ncol = m)
  y[,1] <- xout
  j1 <- floor(n  * (rxo[1] - x[1, 1]) / (x[n, 1] - x[1, 1])) # index of first element
  j1 <- min(n, max(1, j1 - 1))

  j2 <- ceiling(n  * (rxo[2] - x[1, 1]) / (x[n, 1] - x[1, 1]))
  j2 <- max(1, min(n, j2 + 1))
  #cat(j1, j2, "\n")
  if (j1 < j2) {
    for (i in 2:m) {
      y[,i] <- as.vector(approx(x[j1:j2, 1], x[j1:j2, i], xout, ...)$y)
    }
  }
  if (wasdf) y <- as.data.frame(y)
  names(y) <- dimnames(x)[[2]]
  y
} 