#datex <- function (x, how, num = TRUE) {
#  a <- format(x, how)
#  if (num == TRUE) {
#    a <- as.numeric(a)
#  }
#  a
#}

adddoys <- function (dat, col = 1) {
  dat$doy   <- as.numeric(format(dat[, col], "%j"))
  dat$month <- as.numeric(format(dat[, col], "%m"))
  dat$year  <- as.numeric(format(dat[, col], "%Y"))
  dat
}


addyears <- function(ny, uu) {
  uerg <- uu
  for(i in 1:(ny-1)){
    uu2 <- uu
    uu2[,1] <- uu2[,1] + i*365
    uerg <- rbind(uerg, uu2)
  }
  uerg  
}
