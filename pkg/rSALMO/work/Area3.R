library(rSALMO)

data(inputs)
colnames(inputs) <- c("time", "vol", "depth", "dz", "qin", "ased", "srf",
  "iin", "temp", "nin", "pin", "pomin", "zin", "oin", "aver",
  "ad", "au", "eddy", "x1in", "x2in", "x3in", "sf")

nlayers <- 140
  
data(Area3)
## rearrange layers to wide format
alltimes <- unique(inputs[,1])
forcings <- matrix(0, nrow = length(alltimes), ncol = ncol(inputs) * nlayers)
for(i in 1:length(alltimes)){
  forcings[i, ] <- layers2vector(inputs[inputs[ ,1] == alltimes[i], ])
}
## attach the names of forcings as an attribut
attr(forcings, "colnames") <- colnames(inputs)

time <- 100
forc <- approxTime1(forcings, time)

fx <- forc
dim(fx) <- c(22, 140)
fx <- t(fx)
colnames(fx) <- colnames(inputs)

Area3x <- fx[, "vol"] / fx[, "dz"]
Area3x <- c(Area3x, Area3x[length(Area3x)])

