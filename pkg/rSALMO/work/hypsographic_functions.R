data(bautzen_hypso)

hypso <- list() # aggregate hypsographic functions for a lake in one list

hypso$area   <- with(bautzen_hypso, approxfun(level, area))
hypso$volume <- with(bautzen_hypso, approxfun(level, volume))
hypso$level  <- with(bautzen_hypso, approxfun(volume, level))

hypso$vh <- function(level, zmixreal) {
  hypso$volume(level - zmixreal)
}

hypso$ve <- function(level, zmixreal) {
  vh <- hypso$vh(level, zmixreal)
  v  <- hypso$volume(level)
  v - vh
}

hypso$zmix <- function(level, zmixreal) {
  v <- hypso$volume(level)
  vh <- hypso$vh(level, zmixreal)
  a  <- hypso$area(level)
  (v - vh) / a
}

hypso$zhm <- function(level, zmixreal, ah.min = 0.1) {
  vh  <- hypso$vh(level, zmixreal)
  ah  <- hypso$area(level)
  if (ah > ah.min)
    zhm <- vh /ah
  else
    zhm <- 0
  zhm
}
