str(forcings)

cnames <- attr(forcings, "colnames")

ni <- length(cnames)
nlayers <- 140
iTime  <- match("time", cnames)
iDepth <- match("depth", cnames)
iEddy  <- match("eddy",  cnames)
iTemp  <- match("temp",  cnames)

xx <- as.data.frame(matrix(forcings, ncol=ni))
colnames(xx) <- cnames

depth <- forcings[, iDepth + (0:(nlayers - 1) * ni)] 
eddy  <- forcings[, iEddy + (0:(nlayers - 1) * ni)] 
temp <- forcings[, iTemp + (0:(nlayers - 1) * ni)] 
time <- forcings[, iTime + (0:(nlayers - 1) * ni)] 

summary(eddy)

matplot(forcings[,1], depth, type="l")
matplot(forcings[,1], temp, type="l")
matplot(forcings[,1], eddy, type="l")
