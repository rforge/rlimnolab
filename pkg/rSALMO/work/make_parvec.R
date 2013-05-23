#write.table(cc, file="cc.txt", row.names=TRUE, col.names=FALSE, quote=FALSE)

ccx <- read.table("cc.txt")

cc <- as.vector(ccx[,2])
names(cc) <- ccx[,1]
save(cc, file="cc.rda")
     
