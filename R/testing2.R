## Testing from scratch - February 9, 2018

#' Tests to make sure that Fortran calls are working (again)
#' 
#' @details The answers are different due to where sqrt(dt) (?) gets taken into account
#' 
#' @export
#' @useDynLib transfer2
tstWorking <- function(){
  out <- .Fortran("tstWork", x = integer(1))
  out
}

#' Taper Testing
#' 
#' Should get the same answer as from tapering in R
#' 
#' @export
#' @useDynLib transfer2
tstTaperData <- function(){
  d1 <- rnorm(100, mean = 10, sd = 2)
  d2 <- rnorm(200, mean = 5, sd = 0.5)
  
  nw <- 5; k <- 9
  dt1 <- 2; dt2 <- 1; dtRatio <- dt1/dt2
  n1 <- length(d1); n2 <- length(d2)
  nFFT1 <- 2^(floor(log2(n1))+2); nFFT2 <- nFFT1 * dtRatio
  
  out <- .Fortran("tstTaperData", d1 = as.double(d1), d2 = as.double(d2)
                  , n1 = as.integer(n1), n2 = as.integer(n2)
                  , m1 = as.integer(nFFT1), m2 = as.integer(nFFT2)
                  , c1 = complex(nFFT1*k), c2 = complex(nFFT2*k)
                  , dt1 = as.double(dt1), dt2 = as.double(dt2)
                  , k = as.integer(k), nw = as.double(nw))
  
  outtpr1 <- Re(matrix(out$c1, ncol = k, nrow = nFFT1)[1:n1, ])
  outtpr2 <- Re(matrix(out$c2, ncol = k, nrow = nFFT2)[1:n2, ])
  
  v1 <- multitaper::dpss(n = n1, k = k, nw = nw, returnEigenvalues = TRUE)
  v2 <- multitaper::dpss(n = n2, k = k, nw = nw, returnEigenvalues = TRUE)

  tpr1 <- apply(v1$v, 2, "*", d1)
  tpr2 <- apply(v2$v, 2, "*", d2)
  
  print(range(tpr1 - outtpr1))
  print(range(tpr2 - outtpr2))
}

#' Taper Testing
#' 
#' Should get the same answer as from tapering in R
#' 
#' @details Looks good - differences are a result of the implementation of the FFT's between 
#' R and fftpack5 that I'm using in Fortran - differences on the order of 1e-14 at most.
#' 
#' @export
#' @useDynLib transfer2
tstEigenCoef <- function(){
  d1 <- rnorm(100, mean = 10, sd = 2)
  d2 <- rnorm(200, mean = 5, sd = 0.5)
  
  nw <- 5; k <- 9
  dt1 <- 2; dt2 <- 1; dtRatio <- dt1/dt2
  n1 <- length(d1); n2 <- length(d2)
  nFFT1 <- 2^(floor(log2(n1))+2); nFFT2 <- nFFT1 * dtRatio
  
  s1 <- multitaper::spec.mtm(d1, nw = nw, k = k, deltat = dt1, dtUnits = "second", nFFT = nFFT1
                 , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  s2 <- multitaper::spec.mtm(d2, nw = nw, k = k, deltat = dt2, dtUnits = "second", nFFT = nFFT2
                             , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  
  out <- .Fortran("tstEigenCoef", d1 = as.double(d1), d2 = as.double(d2)
                  , n1 = as.integer(n1), n2 = as.integer(n2)
                  , m1 = as.integer(nFFT1), m2 = as.integer(nFFT2)
                  , yk1 = complex( (nFFT1/2+1) * k ), yk2 = complex( (nFFT2/2+1) * k )
                  , dt1 = as.double(dt1), dt2 = as.double(dt2)
                  , k = as.integer(k), nw = as.double(nw))
  
  print("Series1 Diferences: ")
  print("Real:")
  print(paste("Mean =", mean(Re(s1$mtm$eigenCoefs - out$yk1)), " :: "
              , "Median =", median(Re(s1$mtm$eigenCoefs - out$yk1)), " :: "
              , "Range =", paste(range(Re(s1$mtm$eigenCoefs - out$yk1)), collapse = ", ")))
  print(range(Im(s1$mtm$eigenCoefs - out$yk1)))
  
  hist(Re(s1$mtm$eigenCoefs - out$yk1))
  
  
  print("Series2: ")
  print(range(Re(s2$mtm$eigenCoefs - out$yk2)))
  print(range(Im(s2$mtm$eigenCoefs - out$yk2)))
}

#' Adaptive Eigencoefficient Testing
#' 
#' Should get the same answer as from multitaper in R
#' 
#' @details Less good - differences on the order of 1e-7 ... 
#' 
#' @export
#' @useDynLib transfer2
tstWeightedEigenCoef <- function(d1 = NULL, d2 = NULL){
  if (is.null(d1)){
    d1 <- rnorm(100, mean = 10, sd = 2)
  }
  
  if (is.null(d1) || length(d2) != 2*length(d1)){
    d2 <- rnorm(2*length(d1), mean = 5, sd = 0.5)
  }
  
  nw <- 5; k <- 9
  dt1 <- 2; dt2 <- 1; dtRatio <- dt1/dt2
  n1 <- length(d1); n2 <- length(d2)
  nFFT1 <- 2^(floor(log2(n1))+2); nFFT2 <- nFFT1 * dtRatio
  
  s1 <- multitaper::spec.mtm(d1, nw = nw, k = k, deltat = dt1, dtUnits = "second", nFFT = nFFT1
                             , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  s2 <- multitaper::spec.mtm(d2, nw = nw, k = k, deltat = dt2, dtUnits = "second", nFFT = nFFT2
                             , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  
  out <- .Fortran("tstWeightedEigenCoef", d1 = as.double(d1), d2 = as.double(d2)
                  , n1 = as.integer(n1), n2 = as.integer(n2)
                  , m1 = as.integer(nFFT1), m2 = as.integer(nFFT2)
                  , yk1 = complex( (nFFT1/2+1) * k ), yk2 = complex( (nFFT2/2+1) * k )
                  , dt1 = as.double(dt1), dt2 = as.double(dt2)
                  , k = as.integer(k), nw = as.double(nw))
  
  
  wteig1 <- s1$mtm$eigenCoefs * s1$mtm$eigenCoefWt
  wteig2 <- s2$mtm$eigenCoefs * s2$mtm$eigenCoefWt
  
  outyk1 <- out$yk1
  outyk2 <- out$yk2
  
  print("Series 1:")
  print("Real Part:")
  print(paste("Range of spec:", range(Re(wteig1))))
  print(paste("Range = ", paste(range(Re(wteig1 - outyk1)), collapse = ", "), " :: "
              , "Mean = ", mean(Re(wteig1 - outyk1)), " :: "
              , "Median = ", median(Re(wteig1 - outyk1))))
  print(range(Im(wteig1 - outyk1)))
  print("Series 2:")
  print(range(Re(wteig2 - outyk2)))
  print(range(Im(wteig2 - outyk2)))
}

#' @export
#' @useDynLib transfer2
tstAdaptiveWeights <- function(d1 = NULL, d2 = NULL){
  if (is.null(d1)){
    d1 <- rnorm(100, mean = 10, sd = 2)
  }
  
  if (is.null(d1) || length(d2) != 2*length(d1)){
    d2 <- rnorm(2*length(d1), mean = 5, sd = 0.5)
  }
  
  nw <- 5; k <- 9
  dt1 <- 2; dt2 <- 1; dtRatio <- dt1/dt2
  n1 <- length(d1); n2 <- length(d2)
  nFFT1 <- 2^(floor(log2(n1))+2); nFFT2 <- nFFT1 * dtRatio
  v1 <- var(d1); v2 <- var(d2)
  
  
  s1 <- multitaper::spec.mtm(d1, nw = nw, k = k, deltat = dt1, dtUnits = "second", nFFT = nFFT1
                             , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  s2 <- multitaper::spec.mtm(d2, nw = nw, k = k, deltat = dt2, dtUnits = "second", nFFT = nFFT2
                             , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  
  out <- .Fortran("tstAdaptiveWeights", d1 = as.double(d1), d2 = as.double(d2)
                  , n1 = as.integer(n1), n2 = as.integer(n2)
                  , m1 = as.integer(nFFT1), m2 = as.integer(nFFT2)
                  , yk1 = complex( (nFFT1/2+1) * k ), yk2 = complex( (nFFT2/2+1) * k )
                  , dt1 = as.double(dt1), dt2 = as.double(dt2)
                  , k = as.integer(k), nw = as.double(nw)
                  , dk1 = as.double((nFFT1/2+1) * k), dk2 = as.double((nFFT2/2+1) * k)
                  , var1 = as.double(v1), var2 = as.double(v2))
  
  out
}


#' @export
#' @useDynLib transfer2
eigenCoefWt <- function(){
  d1 <- rnorm(100, mean = 10, sd = 2)
  nw <- 5; k <- 9
  dt1 <- 2; dt2 <- 1; dtRatio <- dt1/dt2
  n1 <- length(d1); #n2 <- length(d2)
  nFFT1 <- 2^(floor(log2(n1))+2); #nFFT2 <- nFFT1 * dtRatio
  v1 <- var(d1); #v2 <- var(d2)
  out <- .Fortran("tstEiegenCoefWt", d = as.double(d1)
                  , ndata = as.integer(length(d1))
                  , dt = as.double(dt1)
                  , nw = as.double(nw)
                  , k = as.integer(k)
                  , yk = complex(nFFT1*k)
                  , spec = complex(nFFT1/2+1)
                  , nFFT = as.integer(nFFT1), id = as.integer(1))
}

# weights are exactly the same if I use multitaper yk's, v's, and ev's
# i.e., weight calculation is correct.
#' @export
testWeightsWrapper <- function(d1 = NULL, d2 = NULL){
  if (is.null(d1)){
    d1 <- rnorm(100, mean = 10, sd = 2)
  }
  
  if (is.null(d1) || length(d2) != 2*length(d1)){
    d2 <- rnorm(2*length(d1), mean = 5, sd = 0.5)
  }
  
  nw <- 5; k <- 9
  dt1 <- 2; dt2 <- 1; dtRatio <- dt1/dt2
  n1 <- length(d1); n2 <- length(d2)
  nFFT1 <- 2^(floor(log2(n1))+2); nFFT2 <- nFFT1 * dtRatio
  var1 <- var(d1); var2 <- var(d2)
  
  s1 <- multitaper::spec.mtm(d1, nw = nw, k = k, deltat = dt1, dtUnits = "second", nFFT = nFFT1
                             , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  s2 <- multitaper::spec.mtm(d2, nw = nw, k = k, deltat = dt2, dtUnits = "second", nFFT = nFFT2
                             , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  
  out <- testweights(d1, s1$mtm$eigenCoefs, k = k, dt = dt1, var = var1, nFFT = nFFT1, nw = nw,
                     eval = s1$mtm$dpss$eigen)
  
  dkFort <- matrix(out$dk, ncol = k, nrow = nFFT1/2+1)
  dkMtm <- s1$mtm$eigenCoefWt
  
  print(range(dkMtm - dkFort))
}

#' @export
#' @useDynLib transfer2
tstEigenvals <- function(d1 = NULL, d2 = NULL){
  if (is.null(d1)){
    d1 = rnorm(100, mean = 10, sd = 2)
  }
  
  if (is.null(d1) || length(d2) != 2*length(d1)){
    d2 <- rnorm(2*length(d1), mean = 5, sd = 0.5)
  }
  
  nw <- 5; k <- 9
  dt1 <- 2; dt2 <- 1; dtRatio <- dt1/dt2
  n1 <- length(d1); n2 <- length(d2)
  nFFT1 <- 2^(floor(log2(n1))+2); nFFT2 <- nFFT1 * dtRatio
  var1 <- var(d1); var2 <- var(d2)
  
  s1 <- multitaper::spec.mtm(d1, nw = nw, k = k, deltat = dt1, dtUnits = "second", nFFT = nFFT1
                             , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  s2 <- multitaper::spec.mtm(d2, nw = nw, k = k, deltat = dt2, dtUnits = "second", nFFT = nFFT2
                             , Ftest = TRUE, returnInternals = TRUE, plot = FALSE, centre = "none")
  
  out <- .Fortran("tstEigenvals", n1 = as.integer(n1), n2 = as.integer(n2)
                  , k = as.integer(k)
                  , nw = as.double(nw), m1 = as.integer(nFFT1), m2 = as.integer(nFFT2)
                  , lambda1 = double(k), lambda2 = double(k))
  
  print("Series1:")
  print(s1$mtm$dpss$eigen - out$lambda1, digits = 15)
  print("Series2:")
  print(s2$mtm$dpss$eigen - out$lambda2, digits = 15)
}

#' @export
#' @useDynLib transfer2
tstblockincrement <- function(blockSize = 100, overlap = 0){
  out <- .Fortran("tstblockincrement", block_size = as.integer(blockSize)
                  , overlap = as.double(overlap), block_incr = integer(1))
  
  print(paste("Given overlap:", overlap))
  print(paste("Overlap:", round((blockSize - out$block_incr) / blockSize, 4)))
}

#' @export
#' @useDynLib transfer2
tstcalculatenblocks <- function(ndata = 100, blockSize = 50, blockIncr = 50){
  out <- .Fortran("tstcalculatenblocks", ndata = as.integer(ndata), block_size = as.integer(blockSize)
                  , block_incr = as.integer(blockIncr), nblocks = integer(1))
  
  out$nblocks - length(seq(1, ndata - blockSize + 1, by = blockIncr))
}

#' @export
#' @useDynLib transfer2
tstcalctfwteigen <- function(d1 = NULL, d2 = NULL, npred = 2){
  if (is.null(d1)){
    d1 = rnorm(100, mean = 10, sd = 2)
    d2 <- matrix(c(rnorm(2*length(d1), mean = 5, sd = 0.5)
                   , rnorm(2*length(d1), mean = 33, sd = 0.1)), ncol = 2)
  }
  
  nw <- 5; k <- 9; nw2 = nw
  dt1 <- 2; dt2 <- 1; dtRatio <- dt1/dt2
  n1 <- length(d1); n2 <- length(d2) / npred
  nFFT1 <- 2^(floor(log2(n1))+2); nFFT2 <- nFFT1 * dtRatio
  var1 <- var(d1); var2 <- var(d2)
  nfreqs1 <- nFFT1/2 + 1; nfreqs2 <- nFFT2/2 + 1
  
  out1 <- .Fortran("tstcalctfwteigen", block_incr = as.integer(n1), block_incr2 = as.integer(n2)
                  , block_size = as.integer(n1), block_size2 = as.integer(n2)
                  , nblocks = as.integer(1), d1 = as.double(d1), d2 = as.double(d2)
                  , dt = as.double(dt1), dt2 = as.double(dt2), nw = as.double(nw)
                  , nw2 = as.double(nw2), k = as.integer(k)
                  , nFFT = as.integer(nFFT1), nFFT2 = as.integer(nFFT2)
                  , yk1 = complex(nfreqs1 * k), yk2 = complex(nfreqs2 * k * npred)
                  , ndata = as.integer(n1), ndata2 = as.integer(n2), npred = as.integer(npred))
  yk11 <- matrix(out1$yk1, ncol = k)
  yk21 <- array(out1$yk2, dim = c(nfreqs2, k, npred))
  
  out2 <- .Fortran("tstWeightedEigenCoef", d1 = as.double(d1), d2 = as.double(d2[, 1])
                   , n1 = as.integer(n1), n2 = as.integer(n2)
                   , m1 = as.integer(nFFT1), m2 = as.integer(nFFT2)
                   , yk1 = complex( (nFFT1/2+1) * k ), yk2 = complex( (nFFT2/2+1) * k )
                   , dt1 = as.double(dt1), dt2 = as.double(dt2)
                   , k = as.integer(k), nw = as.double(nw))
  
  yk12 <- matrix(out2$yk1, ncol = k)
  yk22 <- matrix(out2$yk2, ncol = k)
  
  print(range(abs(yk11 - yk12)))
  print(range(abs(yk21[, , 1] - yk22)))
}

#' @export
#' @useDynLib transfer2
tstfindhidx <- function(curCol = 1){
  ind <- array(0, c(16, 16, 3))
  
  ind[4, 6, 1] <- 1
  ind[4, 9, 1] <- 1
  ind[4, 2, 1] <- 1
  ind[9, 1, 1] <- 1
  ind[9, 5, 1] <- 1
  ind[9, 6, 1] <- 1
  ind[9, 7, 1] <- 1
  ind[9, 13, 1] <- 1
  ind[16, 13, 1] <- 1
  ind[16, 11, 1] <- 1
  ind[16, 1, 1] <- 1
  ind[2, 1, 2] <- 1
  ind[5, 5, 2] <- 1
  ind[5, 6, 2] <- 1
  ind[5, 10, 2] <- 1
  ind[1, 10, 3] <- 1
  ind[1, 2, 3] <- 1
  ind[9, 7, 3] <- 1
  ind[10, 7, 3] <- 1
  hIdx <- c()
  totFreqByColByPred <- matrix(0, nrow = 3+1, ncol = 16)
  hPredBreak <- 1
  for (i in 1:3){
    tmp <- (1:16)[apply(ind[, , i], 1, sum) > 0]
    hPredBreak[i+1] <- hPredBreak[i] + length(tmp)
    hIdx[ hPredBreak[i]:(hPredBreak[i+1] - 1) ] <- tmp
    
    totFreqByColByPred[i+1, ] <- apply(ind[, , i], 2, sum)
  }
  
  idxSub <- rep(-1, sum(totFreqByColByPred[, curCol]))
  out <- list()
  for (p in 1:3){
    if (totFreqByColByPred[p+1, curCol] == 0){ next }
    out[[p]] <- .Fortran("tsthfind", col = as.integer(ind[, curCol, p])
                         , hIdx = as.integer(hIdx[hPredBreak[p]:(hPredBreak[p+1]-1)])
                         , ncol = as.integer(16)
                         , nhIdx = as.integer(length(hIdx[hPredBreak[p]:(hPredBreak[p+1]-1)]))
                         , idxSub = as.integer(idxSub[(sum(totFreqByColByPred[1:p, curCol])+1):(sum(totFreqByColByPred[1:(p+1), curCol]))])
                         , nSub = as.integer(sum(totFreqByColByPred[p:(p+1), curCol]))
                         , baseIdx = as.integer(hPredBreak[p]-1))
    idxSub[(sum(totFreqByColByPred[1:p, curCol])+1):(sum(totFreqByColByPred[1:(p+1), curCol]))] <- out[[p]]$idxSub
  }
  
  print(paste("Length of idxSub: ", length(idxSub)))
  print(paste("idxSub = ", paste(idxSub, collapse = ", ")))
  print(hIdx)
  list(out = out, ind = ind, hIdx = hIdx, hPredBreak = hPredBreak, totFreqByColByPred = totFreqByColByPred, idxSub = out$idxSub)
}

#' @export
tsttf <- function(d1 = NULL, d2 = NULL, dt=1, dt2=dt){
  if (is.null(d1)){
    d1 = rnorm(100, mean = 10, sd = 2)
    d2 <- matrix(c(rnorm(2*length(d1), mean = 5, sd = 0.5)
                   , rnorm(2*length(d1), mean = 33, sd = 0.1)), ncol = 2)
    npred <- 2
  }
  
  npred <- dim(d2)[2]
  
  nw <- 5; k <- 9; nw2 = nw
  dt1 <- 2; dt2 <- 1; dtRatio <- dt1/dt2
  n1 <- length(d1); n2 <- length(d2) / npred
  nFFT1 <- 2^(floor(log2(n1))+2); nFFT2 <- nFFT1 * dtRatio
  blockSize <- n1
  
  # d1, d2, ndata = length(d1), ndata2 = dim(d2)[1]
  # , blockSize = ndata, blockSize2 = blockSize, overlap = 0
  # , dt = 1, dt2 = dt, nw = 4, nw2 = NULL, k = 7
  # , nFFT = NULL, nFFT2 = nFFT
  # , freqRange = NULL, maxFreqOffset = 0, nOff = -1
  # , forceZeroFreq = TRUE
  # , sigLevel = 0.99
  # , name1 = "d1", name2 = "d2"
  
  a <- tf(d1 = d1, d2 = as.data.frame(d2), blockSize = length(d1), blockSize2 = 2*length(d1)
          , dt = 2, dt2 = 1, nFFT = nFFT1, nFFT2 = nFFT2)
}