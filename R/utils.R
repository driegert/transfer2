# Determines nw2 such that the effective w1 and w2 are equal.
# This should mean that nw2 is at least as large as nw1
determineNW2 <- function(n1, n2, dt1, dt2, nw1){
  (dt2*n2*nw1) / (dt1*n1)
}