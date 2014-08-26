#' Extract Physical Structure from Input Data
#'
#' The function takes a data-base type table and returns
#'  a list with interpolated matrices for eddy diffusivity,
#'  temperature and depths
#'
#' @param dat    data set with columns ``depth'', ``temp'' and ``eddy''
#' @param times  vector of requested time steps for interpolation
#' @param depths  vector of requested depths for interpolation
#'
#' @return list with matrices ``eddy'', ``temp'' and ``depth''
#'
#' depths <- seq(70, 0, -0.5)
#' times <- unique(dat$time)
#' # example to be continued ...

extract_phys <- function(dat, times, depths) {
  ntimes  <- length(times)
  ndepths <- length(depths)

  tempmat <- eddymat <- depthmat <- matrix(0, ntimes, ndepths)

  for (i in 1:length(times)) {
     time_i <- times[i]
     depth <- with(dat, depth[time == time_i])
     temp  <- with(dat,  temp[time == time_i])
     eddy  <- with(dat,  eddy[time == time_i])
     tempmat[i,] <- approx(depth, temp, depths, rule = 2)$y
     eddymat[i,] <- approx(depth, eddy, depths, rule = 2)$y
     depthmat[i,] <- depths
  }
  list(eddy = eddymat, temp = tempmat, depth = depthmat)
}

