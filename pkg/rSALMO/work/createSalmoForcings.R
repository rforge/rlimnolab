## code fragment, just at the very beginning

createSalmoForcings <- function() {

  forcings <- data.frame(
    time   = 0:365,
    vol    = vol,
    depth  = depth,
    dz     = dz,
    qin    = qin,
    ased   = ased,
    srf    = 0,
    iin    = iin,
    temp   = temp,
    nin    = nin,
    pin    = pin,
    pomin  = pomin,
    zin    = zin,
    oin    = oin,
    aver   = aver,
    ad     = 0,
    au     = 0,
    diff   = diff,
    x1in   = X1in,
    x2in   = X2in,
    x3in   = X3in,
  )
  forcings
}
