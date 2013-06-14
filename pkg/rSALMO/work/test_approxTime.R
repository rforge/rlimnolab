x <- data.frame(time = c(1,2,4,5,7,11),
              y1 = runif(6),
              y2 = runif(6)
)

x$time <- x$time * 0.3

matplot(x[,1], x[,-1], type="l")

#ret <- approxTime1(x, 7.5)
#matpoints(ret[1], t(ret[-1]), pch=16)


(ret <- approxTimeEq(x, 1.6))
matpoints(ret[1], ret[-1], pch=16)


size <- 1e2
x <-  x <- data.frame(time = (1:size) * 1.5 - 50,
              y1 = runif(size),
              y2 = runif(size)
)

matplot(x[,1], x[,-1], type="l")
for (i in 1:10) {
  ret <- approxTimeEq(x, xo <- runif(3, min=min(x$time), max=max(x$time)))
  print(ret)
  matpoints(ret[,1], ret[,-1], pch=16)
}

system.time(
  for (i in 1:10)
    ret <- approxTimeEq(x, xo <- runif(1, min=min(x$time), max=max(x$time)))
)

system.time(
  for (i in 1:10)
    ret <- approxTime(x, -1)
)


