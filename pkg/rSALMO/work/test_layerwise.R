
nstates <- 10
nlayers <- 3

swise <- rep(1:nstates, each=nlayers) + rep((1:3)/10, nstates)

matrix(swise, nrow=nlayers, ncol=nstates)


lwise <- as.vector(t(matrix(swise, nrow=nlayers, ncol=nstates))) # layer wise
lwise

swise2 <- as.vector(t(matrix(lwise, nrow=nstates, ncol=nlayers))) #state wise
swise2

swise == swise2
