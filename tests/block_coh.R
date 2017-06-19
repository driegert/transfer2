library(fields)

bou <- readRDS(file="~/school_lab/phd_proposal/assets/data/geomag/boulder_1999-2012_15min_final.rds")
for (i in 0:15){
  numb <- i
  n2 <- 1
  idxUse <- (numb*28800 + 1):(numb*28800 + n2*28800)
  d1 <- bou$Y[idxUse] - lm(bou$Y[idxUse] ~ idxUse)$fitted.value
  d2 <- bou$Z[idxUse] - lm(bou$Z[idxUse] ~ idxUse)$fitted.value
  
  startTime <- proc.time()
  a <- coherence(d1, d2, ndata = length(d1), blockSize = 28800, overlap = 0, 
                 , dt = 900, nw = 4, k = 7
                 , freqRange = c(60e-6, 120e-6), maxFreqOffset = 30e-6, calcType = 1)
  proc.time() - startTime
  range(Re(a$coh))
  
  # coh.offset <- apply(abs(a$coh)^2, 1, mean)
  coh.offset <- apply(Re(a$coh), 1, mean)
  
  cpd <- 1/(24*3600)
  mcpd <- c(rev(seq(0, -3*cpd, by = -cpd/2)), seq(0, 3*cpd, by = cpd/2))
  par(mar = c(4,4,1,1))
  plot(a$offset*1e6, coh.offset, type='l')
  abline(v = 1e6*mcpd, lty = 3, col = 'orange')
}
# par(mar = c(4,4,1,1))
# image.plot(x = a$bandfreq, y = a$offset, z = t(a$coh))

par(mar = c(4,4,1,1))
plot(d1[1:28800], type='l')
plot(d2[1:28800], type='l')



dt = 900
freqRange = c(60, 120)*1e-6
b <- coherence(d1, d2, ndata = length(d1), blockSize = 28800, overlap = 0, 
               , dt = dt / dt, nw = 4, k = 7
               , freqRange = freqRange * dt, maxFreqOffset = dt * 30e-6)
b.offset <- apply(abs(b$coh)^2, 1, mean)
