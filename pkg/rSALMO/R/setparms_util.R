## ===============================================================
## handling of model parameters
## local helper functions
## ===============================================================

list2matrix <- function(x) {
  if (!is.list(x)) stop("Argument is not a list.")
  
  ncol <- unique(sapply(x, function(element) length(element)))
  if (length(ncol) > 1) 
    stop("Data structure does not have unique element length.")
  ret <- matrix(unlist(x), ncol = ncol, byrow = TRUE)
  nam <- names(x)
  row.names(ret) <- nam
  
  # convert to vector if matrix has only one column
  if(ncol(ret) == 1) 
    ret <- ret[,1]
  
  ret
}

## change model parameters in a named vector or matrix
setparms_matrix <- function() {
  if (is.matrix(parms)) {
    xnames <- rownames(parms)
  } else if (is.vector(parms)) {
    xnames <- names(parms)
  } else {
    stop("Don't know how to handle obj of class ", class(parms))
  }
  
  ## modify vector or matrix component
  xnames   <- toupper(xnames)
  ynames   <- toupper(pnames)
  
  matching <- ynames %in% xnames
  
  missing   <- pnames[!(matching)]
  if (length(missing > 0))
    warning("Parameter/s) ", paste(missing, collapse =", "), " not found")
  
  if(is.matrix(parms)) {
    parms[na.omit(pmatch(ynames, xnames)), col] <- values[which(matching)]  
  } else {
    parms[na.omit(pmatch(ynames, xnames))] <- values[which(matching)]  
  }
  parms
}

## change model parameters in a list of vectors and matrices
setparms_list <- function(old, new) {
  modify_slot <- function(slot) {
    ## columns of a matrix
    ncol <- dim(new[[slot]])[2]
    nam <- row.names(new[[slot]])
    ## or 1 if new[[slot]] is a vector
    if (is.null(ncol)) {
      ncol <- 1
      nam <- names(new[[slot]])
    }
    x <- setparms(old[[slot]], nam, new[[slot]], 1:ncol)
    x
  }
  new <- lapply(new, list2matrix)
  
  both    <- list(old, new)
  matched <- unique(unlist(lapply(both, names)))
  names(matched) <- matched
  
  ret <- old
  ret <- lapply(matched, modify_slot)
  ret
}

## main function for changing model parameters
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
      parms[[slot]] <- setparms(parms[[slot]], pnames = pnames, values = values, col = col, slot = NA)
    } else {
      warning("Slot not found")  
    }
    ## if 1st argument is a vector or matrix  
  } else {
    parms <- setparms_matrix(parms, values)
  }
  parms
}

