#' Reformat Layer Strcuture from Long Format to Wide Format
#'
#' The function takes a data-base type table and returns
#'  a table in wide format.
#'
#'
#' @param x matrix or data frame that is to be reformatted
#' @param rowname name of the variable that is used as the row name
#' @param colname name of the variable that is used as the column name
#'
#' @return matrix
#'
#' \dontrun{
#' forcings <- read.table("forcings.txt", header=TRUE)
#' x <- as.matrix(forcings)
#' xx <- layers_long2wide(x, "time", "depth")
#' } % end dontrun



layers_long2wide <- function(x, rowname, colname) {
  if (is.data.frame(x)) x <- as.matrix(x)
  columns <- unique(x[, colname])
  rows    <- unique(x[, rowname])
  ncols   <- length(columns)
  nrows   <- length(rows)
  res <- matrix(0, nrow = nrows, ncol = ncol(x) * ncols)
  for(i in 1:nrows)
    res[i, ] <- layers2vector(x[x[ ,rowname] == rows[i], ])

  return(res)
}

