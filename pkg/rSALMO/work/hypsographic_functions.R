data(bautzen_hypso)


set_hypso_functions <- function(table) {
  level  <- with(table, approxfun(volume, level))
  area   <- with(table, approxfun(level, area))
  volume <- with(table, approxfun(level, volume))
  
  vh <- function(level, zmixreal) {
    volume(level - zmixreal)
  }
  
  ve <- function(level, zmixreal) {
    vh <- vh(level, zmixreal)
    v  <- volume(level)
    v - vh
  }
  
  zmix <- function(level, zmixreal) {
    v <- volume(level)
    vh <- vh(level, zmixreal)
    a  <- area(level)
    (v - vh) / a
  }
  
  zhm <- function(level, zmixreal, ah.min = 0.1) {
    vh  <- vh(level, zmixreal)
    ah  <- area(level)
    if (ah > ah.min)
      zhm <- vh /ah
    else
      zhm <- 0
    zhm
  }
  
  list(level = level, area = area, volume = volume,
       ve = ve, vh = vh, zmix = zmix, zhm = zhm)
}

hypso <- set_hypso_functions(bautzen_hypso)

