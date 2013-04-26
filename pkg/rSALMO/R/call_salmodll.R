#' Call SALMO Shared Library
#'
#' Call the shared library of SALMO
#'
#'
#' @param cfunc  string, name of the C function to be called
#' @param nOfVar vector with number of variables
#' @param cc     vector of constants for the model
#' @param pp     matrix of phytoplankton parameters
#' @param uu     input vector (environmental conditions)
#' @param xx     state vector
#'
#' @return list, contains past state and derivatives separated into
#'               source and sink term

call_salmodll <- function(cfunc, nOfVar, cc, pp, uu, xx) {

  ## allocate memory
  dxq <- dxs <- numeric(length(xx))
  ## call SALMO core function
  ret  <- .C(cfunc, as.integer(nOfVar), cc = as.double(cc),
             pp = as.double(pp), uu = as.double(uu), xx = as.double(xx),
             dxq = as.double(dxq), dxs = as.double(dxs))
  list(ret$xx, ret$dxq, ret$dxs)
}
