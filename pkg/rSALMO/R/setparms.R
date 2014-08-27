#' Set Model Parameters
#' 
#' Robustly replace model parameters by ignoring capitalisation and order of
#' replacment.
#' 
#' 
#' @param parms data object containing model parameters as a vector, a matrix
#'   or a non-nested list of vectors and matrices
#' @param values vector or list of replacement values
#' @param pnames if \code{values} is a vector or matrix: 
#'   vector of parameter names that will be replaced
#' @param col column in matrix-like parameter objects
#' @parms slot character, name of a single list element of \code{parms}
#' @return the manipulated \code{obj} (vector, matrix or list)
#' @examples
#' 
#' data(cc)
#' cc <- setparms(cc, c("YZP","dummy", "EPSMIN", "dtt"), c(0.7, 66, 0.77, 0.1))
#' cc
#'
#' data(pp)
#' pp <- setparms(pp, c("TOPTX", "vs"), c(28, 0.22), c(1, 3))
#' # pp <- setparms(pp, "KP", c(1.4, 1.5, 1.6), 1:3) # not yet possible
#' 
#' @rdname setparms
#' @export setparms
#'
setparms <- function(parms, values, pnames = NULL, col = 1, slot = NA) {
  #if (length(pnames) != length(values)) 
  #  stop("Length of pnames and values do not match.")
  parvec <- parms # will be extended to other objects
  
  
  if (is.list(parms)) {
    ## second argument is a list of changed parameters
    if (is.list(values)) {
      parms <- setparms_list(parms, values)
      ## if only first argument is a list, call setparms recursively one level
    } else if (!is.null(parms[[slot]])) {
      parms[[slot]] <- setparms(parms[[slot]], values = values, pnames = pnames, col = col, slot = NA)
    } else {
      cat(slot, "\n")
      warning("Slot not found")  
    }
    ## if 1st argument is a vector or matrix  
  } else {
    parms <- setparms_matrix(parms, values, pnames, col = col)
  }
  parms
}


## old version -----------------------------------------------------------------
# setparms <- function(obj, pnames, values, col = 1, slot = NA) {
#   if (length(pnames) != length(values)) 
#     stop("Length of pnames and values do not match.")
#   parvec <- obj # will be extended to other objects
#   
#   ## if list, call setparms recursively one level
#   if (is.list(obj)) { 
#     if (!is.null(obj[[slot]])) {
#       obj[[slot]] <- setparms(obj[[slot]], pnames = pnames, values = values, col = col, slot = NA)
#     } else {
#       warning("Slot not found")  
#     }
#   } else {
#     if (is.matrix(obj)) {
#       xnames <- rownames(obj)
#     } else if (is.vector(obj)) {
#       xnames <- names(obj)
#     } else {
#       stop("Don't know how to handle obj of class ", class(obj))
#     }
# 
#     ## modify vector or matrix component
#     xnames   <- toupper(xnames)
#     ynames   <- toupper(pnames)
#     
#     matching <- ynames %in% xnames
#     
#     missing   <- pnames[!(matching)]
#     if (length(missing > 0))
#       warning("Parameter/s) ", paste(missing, collapse =", "), " not found")
#     
#     if(is.matrix(obj)) {
#       obj[na.omit(pmatch(ynames, xnames)), col] <- values[which(matching)]  
#     } else {
#       obj[na.omit(pmatch(ynames, xnames))] <- values[which(matching)]  
#     }
#   }
#   obj
# }
