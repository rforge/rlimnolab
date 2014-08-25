#' Convert Layer Structure to Vector
#' 
#' These functions are used to reorganize the state vector state-wise resp.
#' layer-wise.
#' 
#' Functions \code{arrangeStateWise} and \code{arrangeLayerWise} rearrange
#' state vectors from layer wise (ABCABCABC) to state wise (AAABBBCCC) and vice
#' versa, function \code{layers2vector} creates a state vector from a matrix
#' (layers x variables).
#' 
#' @aliases arrangeLayerWise arrangeStateWise layers2vector
#' @param x state vector resp. state matrix
#' @param nstates number of state variables
#' @param nlayers number of layers
#' @return rearranged vector of state variables
#'
#' @rdname layers2vector
#' @export arrangeStateWise
arrangeStateWise <- function(x, nstates, nlayers) {
  ndx <- rep(1:nstates, each = nlayers) + rep((0:(nlayers - 1)), nstates) * nstates
  x[ndx]
}

#' @rdname layers2vector
#' @export arrangeLayerWise
arrangeLayerWise <- function(x, nstates, nlayers) {
  arrangeStateWise(x, nlayers, nstates)
}

#' @rdname layers2vector
#' @export layers2vector
layers2vector <- function(x) {
  nlayers <- nrow(x)
  nstates <- ncol(x)
  arrangeLayerWise(as.vector(x), nstates, nlayers)
}
