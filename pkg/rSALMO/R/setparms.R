#' Set Model Parameters
#' 
#' Robustly replace model parameters by ignoring capitalisation and order of
#' replacment.
#' 
#' 
#' @param obj data object containing model parameters (currently only vector
#' and matrix implemented)
#' @param pnames vector of parameter names that shoudl be replaced
#' @param values vector of replacement values
#' @param col column in matrix-like parameter objects
#' @parms slot character id of a list element if \code{obj} is a list
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
#' @export setparms

setparms <- function(obj, pnames, values, col = 1, slot = NA) {
  if (length(pnames) != length(values)) 
    stop("Length of pnames and values do not match.")
  parvec <- obj # will be extended to other objects
  
  ## if list, call setparms recursively one level
  if (is.list(obj)) { 
    if (!is.null(obj[[slot]])) {
      obj[[slot]] <- setparms(obj[[slot]], pnames = pnames, values = values, col = col, slot = NA)
    } else {
      warning("Slot not found")  
    }
  } else {
    if (is.matrix(obj)) {
      xnames <- rownames(obj)
    } else if (is.vector(obj)) {
      xnames <- names(obj)
    } else {
      stop("Don't know how to handle obj of class ", class(obj))
    }

    ## modify vector or matrix component
    xnames   <- toupper(xnames)
    ynames   <- toupper(pnames)
    
    matching <- ynames %in% xnames
    
    missing   <- pnames[!(matching)]
    if (length(missing > 0))
      warning("Parameter/s) ", paste(missing, collapse =", "), " not found")
    
    if(is.matrix(obj)) {
      obj[na.omit(pmatch(ynames, xnames)), col] <- values[which(matching)]  
    } else {
      obj[na.omit(pmatch(ynames, xnames))] <- values[which(matching)]  
    }
  }
  obj
}
