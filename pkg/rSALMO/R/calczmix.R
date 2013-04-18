calczmix <- function(temp, depths) {
  nlayers    <- length(temp)
  t_tresh    <- 1 * (depths[2] - depths[1])
  idzmix     <- which(abs(diff(temp)) >= t_tresh)[1]
  if(is.na(idzmix)) idzmix <- nlayers
  zres       <- depths[idzmix]
  ## return values
  list(idzmix = idzmix, zres = zres)
}
