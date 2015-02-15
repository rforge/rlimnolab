#' Make Input Vector from List of Input Data
#'
#' Interpolates consistent vector of input data for a given time point from 
#' a list that contains data of different spatial and temporal resulution.
#' 
#' @aliases makeInputVector
#'
#' @param inputs list of input data; must contain vectors depth and time
#' @param x vector of requested time steps for interpolation
#'
#'
#' @return vector of interpolated data in the order abcabcabc.
#'
#' @export makeInputVector
#' 
makeInputVector <- function(inputs, x) {

  time <- inputs$time
  i  <- findIndexEq(time, x, rule=2)
  x1 <- max(1, time[i - 1])
  x2 <- time[i]
  ntime  <- length(time)
  ndepth <- ncol(inputs$depth)

  ## todo: make these workarounds obsolete
  inputs$iin <- rep(inputs$iin, each=12)
  inputs$vmatsedi <- NULL
  ## end workaround
  
  ## int2: helper function that is not exported
  int2 <- function(i, xout, x1, x2, y) {
    if ((xout < x1) | (x2 < xout)) {
      y <- y[i,]
    } else {
      y1 <- y[i - 1, ]
      y2 <- y[i, ]
      y <- y1 + (y2 - y1) * (xout - x1)/(x2 - x1)  
    }
    y
  }
  
  
  interpolateInput <- function(y) {
    ## identify structure of data element
    ret <- NULL
    if(is.vector(y)) {
      nr <- length(y)
      nc <- 1
    } else {
      nr <- nrow(y)
      nc <- ncol(y)
    }
    ## interpolate resp. recycle data element
    if (nc == 1) {
      if (nr == 1) {                                   # scalar
        ret <- rep(y, ndepth)
      } else if (nr == ndepth) {                       # depth vector
        ret <- y 
      } else if (nr == ntime) {                        # time dependent
        y <- as.matrix(y, ncol=1)
        ret <- rep(int2(i, x, x1, x2, y), ndepth)
      }
    }  else if ((nr == ntime) & (nc == ndepth)) {      # time x depth matrix
      ret <- int2(i, x, x1, x2, y)
    } else {                                           # otherwise
      ret <- NULL# omit
      
    }
    ret
  }
  abc <- sapply(inputs, interpolateInput, USE.NAMES=FALSE)
  abcabc <- as.vector(t(abc))
  #aabbcc <- as.vector(abc)
  #list(abc=abc, abcabc=abcabc, aabbcc=aabbcc)
  abcabc
}


## ToDo / ideas:
## - makeInputVector a non-local function
## - rename arguments
## - make it a closure (or object method) that classifies inputs according 
##   to a time x space resolution type
## - allow different spatial and temporal resolutions
