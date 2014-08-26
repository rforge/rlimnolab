library(rSALMO)

testdat <- read.table("../xlocal/results.dat", header=TRUE)

# !!!
# rSALMO inputs sind stundenweise aufgelöst
# Outputs des physikalischen Modells liegen stündlich vor
# ==> Performance-Probleme bei Interpolation ?

### Algorithmus zum kopieren von testdat$nue nach inputs$eddy ???

# - Datum ermitteln
# - jeweils kleinste Tiefe ist Oberfläche --> schwierig
# - deshalb von unten nach oben interpolieren

## Methode: "abtasten"; Problem systematischer Fehler da Tag/Nacht

tx <- round(testdat$stunde/24, 2)
txu <- unique(tx)

data(inputs)
colnames(inputs) <- c("time", "vol", "depth", "dz", "qin", "ased", "srf", 
                      "iin", "temp", "nin", "pin", "pomin", "zin", "oin", "aver", 
                      "ad", "au", "eddy", "x1in", "x2in", "x3in", "sf")

#inputs  <- as.data.frame(inputs)
#itime   <- inputs$time
itime <- inputs[, "time"]
uitime <- unique(itime)

for (ty in uitime) {
  cat(ty, "\n")
  ## ToDo: Mittelwertbildung über den Tag
  txx <- txu[which.min(abs(30 - txu))] # nearest time
  
  tmp <- testdat[tx == txx, ]
  xdepth <- tmp$tiefe
  xeddy  <- tmp$nue
  xtemp  <- tmp$temp
  
  ydepth <- inputs[itime == ty, "depth"]

  # Layer von Oberflaeche her gezaehlt; Vorsicht weg. bottom shear!!
  #inputs[itime == ty, "eddy"] <- approx(xdepth - min(xdepth), xeddy, ydepth, rule=2)$y
  #inputs[itime == ty, "temp"] <- approx(xdepth - min(xdepth), xtemp, ydepth, rule=2)$y
  xx <- cbind(xdepth - min(xdepth), xeddy, xtemp)
  
  inputs[itime == ty, c("eddy", "temp")] <- approxTimeEq(xx, ydepth, rule=2)[,-1]
}

#######################################################################################

infl <- read.table("../xlocal/ste_inflows.dat", header=TRUE)
#str(infl)

iin <- infl$iin

cnam <- names(infl)

gw1 <- infl[c(1, grep(".gw1", cnam))]
rai <- infl[c(1, grep(".rai", cnam))]
npp <- infl[c(1, grep(".npp", cnam))]
npp <- npp[-11] # remove qout

names(gw1) <- gsub(".gw1", "", names(gw1))
names(rai) <- gsub(".rai", "", names(rai))
names(npp) <- gsub(".npp", "", names(npp))

gw1 <- cbind(gw1, inflow = "gw1")
rai <- cbind(rai, inflow = "rai")
npp <- cbind(npp, inflow = "npp")

rm(rai, gw1, npp, infl)

inflow <- rbind(gw1, rai, npp)
inflow <- inflow[inflow$time %in% unique(inputs[,1]),] # nur 3 Jahres-Zeitabschnitt







