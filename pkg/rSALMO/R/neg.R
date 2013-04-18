neg <-
function(x, nstates=140) {
  a <- x / nstates
  b <- x %% nstates
  print(a)
  print(b)
  list(a, b)
}
