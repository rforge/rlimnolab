x <- data.frame(time = c(1,2,4,5,7,11),
              y1 = runif(6),
              y2 = runif(6)
)


matplot(x[,1], x[,-1], type="l", xlim=c(-2,15))

approxTime1(x, 14)

x
approxTime1(x, 2.5, rule=2)