voltss <-
function (z, sp = 438.8) {
  Vzmax <- 46.1833
  zz <- z + 438.8 - sp
  Vo <- 10^6 * 22.364
  Vo * (1 - zz/Vzmax)^3
}
