#' Plot Inputs of a 1D Model
#'
#' 
#'
#' @param x input matrix (forcing data)
#' @param variable variable name to be plotted
#' @param plot boolean whether the plot shoudl be shown \code{TRUE}
#' @param \dots further arguments passed to image2D
#' @return subset containing only the requested variable, useful for further analysis
#'
#' #example
#  # x <- plot_inputs(forcings, "eddy")
#  # image2D(log(x[, ncol(x):1]))
#'
#' @export plot_inputs
#'

plot_inputs <- function(x, variable, plot = TRUE, ...) {
  nam <- attr(x, "colnames")
  ndx <- which(nam == variable)
  ndx <- seq(ndx, ncol(x), length(nam))
  nr <- nrow(x)
  x <- x[,ndx]
  if (plot) image2D(x, ...)
  invisible(x)
}
