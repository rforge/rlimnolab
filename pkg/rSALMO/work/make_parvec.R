#write.table(cc, file="cc.txt", row.names=TRUE, col.names=FALSE, quote=FALSE)

ccx <- read.table("cc.txt")

## convert data frame to named vector

## original method
cc <- as.vector(ccx[,2])
names(cc) <- ccx[,1]

## simplified, but not yet tested
#cc <- structure(as.vector(ccx[,2]), names = as.character(ccx[,1]))

save(cc, file="cc.rda")

