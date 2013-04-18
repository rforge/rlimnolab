header <- c("t", "s", "v", "ve", "vh", "zmixreal", "zmix", "zhm", "qin", "qhin", "qout", "qhout", "srf", "iin", "temp", "temph", "pin", "phin", "nin", "nhin", "pomin", "pomhin", "xin1", "xhin1", "xin2", "xhin2", "zin", "zhin", "oin", "ohin", "ae", "ah")

dat <- read.table("uref-1997", header=FALSE)

names(dat) <- header

write.table(dat, "bautzen1997.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)