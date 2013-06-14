## special versions of approx:
##  approxTime:  interpolation of complete rows of a matrix or data frame
##  approxTime1: special case with one row only (slightly faster)
##  approxTimeEq: for equidistant time
##  ToDo: merge these functions and make closures

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
