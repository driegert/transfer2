#' Estimates the transfer functions
#' 
#' Estimates a transfer function using frequency offsets between a response and multiple inputs
#' 
#' @param d1 A real vector containing the response.
#' @param d2 A \code{data.frame} whose columns contain the predictors.
tf <- function(d1, d2, ndata = length(d1), ndata2 = length(d2)
               , blockSize = ndata, blockSize2 = blockSize, overlap = 0
               , dt = 1, dt2 = dt, nw = 4, nw2 = NULL, k = 7
               , nFFT = NULL, nFFT2 = nFFT
               , freqRange = NULL, maxFreqOffset = 0, nOff = -1
               , sigLevel = 0.99
               , name1 = "d1", name2 = "d2"){
  # You need to be hella careful with nFFT and nFFT2:
  # nFFT2 / nFFT needs to be an integer, otherwise... things will go poorly
  # One (me) assumes that you will be using powers of 2 ... 
  if (dt < dt2){
    stop("d1 is the central frequency series.  You should make the faster sampled series d2.")
  }
  
  if (is.null(nFFT) || nFFT < blockSize) {
    nFFT <- 2^(floor(log2(blockSize))+2)
  }
  
  if (is.null(nFFT2) || nFFT2 < blockSize2){
    nFFT2 <- 2^(floor(log2(blockSize2))+2)
  }
  
  fRatio <- (nFFT2 * dt2) / (nFFT * dt)
  
  if (fRatio - trunc(fRatio) != 0){
    stop("You should ensure that your choices of nFFT and nFFT2 work with dt and dt2.\n 
         i.e., df2 / df1 should be an integer.")
  }
  
  # if nw2 is null, make the effective W's the same
  if (is.null(nw2)){
    w1 <- nw / blockSize
    nw2 <- w1*blockSize2
  }
  
  df <- 1 / (dt*nFFT)
  
  # this if-statement is here because there were issues when freqRange == NULL
  if (is.null(freqRange)){
    warnings("freqRange == NULL: Setting freqRange to positive band and maxFreqOffset to 0.")
    freqRange <- c(0, 1/(2*dt))
    maxFreqOffset <- 0
    
    freqRangeIdx <- c(1, nFFT/2+1)
  } else {
    freqRangeIdx <- c(max(1, floor(freqRange[1] / df)), ceiling(freqRange[2]/df))
  }
  
  # print(paste0("ndata: ", ndata))
  maxOffIdx <- ceiling(maxFreqOffset / df)
  freq <- seq(0, 1/(2*dt), 1/(nFFT*dt))
  
  numRow = 2*maxOffIdx + 1
  numCol = freqRangeIdx[2] - freqRangeIdx[1] + 1
  nblocks = ndata / (blockSize * (1 - overlap))
  
  # determine which frequencies get used
  level <- qnorm(sigLevel) / sqrt(nblocks)
  
  # calculate the pairwise MSC's
  ind <- array(data = 0, dim = c(numRow, numCol, ncol(d2)))
  for (i in 1:ncol(d2)){
    out <- .Fortran("callblockcoh", d1 = as.double(d1), d2 = as.double(d2[, i])
                    , ndata = as.integer(ndata), ndata2 = as.integer(ndata2)
                    , block_size = as.integer(blockSize), block_size2 = as.integer(blockSize2)
                    , overlap = as.double(overlap), dt = as.double(dt), dt2 = as.double(dt2)
                    , nw = as.double(nw), nw2 = as.double(nw2), k = as.integer(k)
                    , nFFT = as.integer(nFFT), nFFT2 = as.integer(nFFT2)
                    , fRatio = as.integer(fRatio)
                    , coh = complex(numRow*numCol), cohnrow = as.integer(numRow)
                    , cohncol = as.integer(numCol)
                    , freq = as.double(freq)
                    , offsets = double(numRow), freq_range_idx = as.integer(freqRangeIdx)
                    , max_freq_offset_idx = as.integer(maxOffIdx)
                    , calc_type = as.integer(1)
                    , conv_msc2norm = as.integer(1)
                    , is_forward = as.integer(1))
    out2 <- .Fortran("callMscIndicator", msc = as.double(out$coh), nrow = as.integer(numRow)
                     , ncol = as.integer(numCol), ind = integer(numRow * numCol)
                     , level = as.double(level), nOff = as.integer(nOff))
    
    ind[, , i] <- matrix(out2$ind, nrow = numRow, ncol = numCol)
  }
  
  ind[ind > 0] <- 1 # just replace this to get number of frequencies used per central freq.
  
  ## need to take into account the use of the 0 frequency (optionally included?)
  hIdx <- c()
  hPredBreak <- 1
  for (i in 1:ncol(d2)){
    tmp <- (1:numRow)[apply(ind[, , i], 1, sum) > 0]
    hPredBreak[i+1] <- hPredBreak[i] + length(tmp) - 1
    hIdx[ hPredBreak[i]:hPredBreak[i+1] ] <- tmp
  }
  
  totFreqByCol <- apply(2, hIdx, sum)
  
  # hFreq <- out$offsets[hIdx]
  
  # hOut <- .Fortran("calltf", )
  
}