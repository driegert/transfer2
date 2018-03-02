#' Estimates the transfer functions
#' 
#' Estimates a transfer function using frequency offsets between a response and multiple inputs
#' 
#' @param d1 A real vector containing the response.
#' @param d2 A \code{data.frame} whose columns contain the predictors.
#' 
#' @export
#' @useDynLib transfer2
tf <- function(d1, d2, ndata = length(d1), ndata2 = dim(d2)[1]
               , blockSize = ndata, blockSize2 = blockSize, overlap = 0
               , dt = 1, dt2 = dt, nw = 4, nw2 = NULL, k = 7
               , nFFT = NULL, nFFT2 = nFFT
               , freqRange = NULL, maxFreqOffset = 0, nOff = -1
               , forceZeroFreq = TRUE
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
  
  # M1 / (1/dt1) == a * M2 / (1/dt2), solve for 'a'
  fRatio <- (nFFT2 * dt2) / (nFFT * dt)
  
  if (fRatio - trunc(fRatio) != 0){
    stop("You should ensure that your choices of nFFT and nFFT2 work with dt and dt2.\n 
         i.e., df2 / df1 should be an integer.")
  }
  
  # if nw2 is null, make the effective W's the same
  # w1 <- nw / blockSize
  if (is.null(nw2)){
    # nw2 <- w1*blockSize2
    nw2 <- determineNW2(n1 = blockSize, n2 = blockSize2, dt1 = dt, dt2 = dt2, nw1 = nw)
  }
  
  df <- 1 / (dt*nFFT)
  
  # this if-statement is here because there were issues when freqRange == NULL
  if (is.null(freqRange)){
    warnings("freqRange == NULL: Setting freqRange to full positive band and maxFreqOffset to 0.")
    freqRange <- c(0, 1/(2*dt))
    maxFreqOffset <- 0
    
    freqRangeIdx <- c(1, nFFT/2+1)
  } else {
    freqRangeIdx <- c(max(1, floor(freqRange[1] / df)), ceiling(freqRange[2]/df))
  }
  
  # print(paste0("ndata: ", ndata))
  maxOffIdx <- floor(maxFreqOffset / df)
  freq <- seq(0, 1/(2*dt), 1/(nFFT*dt))
  
  numRow <- 2*maxOffIdx + 1
  numCol = freqRangeIdx[2] - freqRangeIdx[1] + 1
  zeroFreq <- maxOffIdx + 1
  zeroDeadZone <- ceiling((w1 / 3) / df)
  nblocks <- ndata / (blockSize * (1 - overlap))
  
  # determine which frequencies get used
  level <- qnorm(sigLevel) / sqrt(nblocks)
  
  if (forceZeroFreq & (maxOffIdx < zeroDeadZone)){
    stop("maxFreqOffset should be larger than (B1)/3 in order to use the zero offset frequency.")
  }
  
  # Just need to do a "normal" transfer function regression - write some code in Fortran.
  if (maxOffIdx == 0){
    
  }
  
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
    if (maxFreqOffset == 0)
    out2 <- .Fortran("callMscIndicator", msc = as.double(Re(out$coh)), nrow = as.integer(numRow)
                     , ncol = as.integer(numCol), ind = integer(numRow * numCol)
                     , level = as.double(level), nOff = as.integer(nOff))
    
    ind[, , i] <- matrix(out2$ind, nrow = numRow, ncol = numCol)
    
    ## something about 
    
    # make the zero-offset a "significant" frequency and implement a deadzone around it.
    if (forceZeroFreq & (maxOffIdx > zeroDeadZone)){
      ind[(zeroFreq - zeroDeadZone):(zeroFreq+zeroDeadZone), , i] <- 0
      ind[zeroFreq, , i] <- 1
    }
  }
  
  # rather than a ranking, just replace with 1's as indicators
  ind[ind > 0] <- 1 # just replace this to get number of frequencies used per central freq.
  
  # this whole piece gives the significant frequencies for all predictors in a single
  # vector, hIdx
  # hPredBreak gives the start and end+1 for where each predictor series lives
  # e.g., if pred 1 had 4, 9, 16; pred 2 had 2, 5; and pred 3 had 1, 9, 10
  # hIdx = c(1,4,9,16, 2,5, 1,9,10); hPredBreak = c(1, 4, 6, 9)
  hIdx <- c()
  totFreqByColByPred <- matrix(0, nrow = ncol(d2)+1, ncol = numCol)
  hPredBreak <- 1
  for (i in 1:ncol(d2)){
    tmp <- (1:numRow)[apply(ind[, , i], 1, sum) > 0]
    hPredBreak[i+1] <- hPredBreak[i] + length(tmp)
    hIdx[ hPredBreak[i]:(hPredBreak[i+1] - 1) ] <- tmp
    
    totFreqByColByPred[i+1, ] <- apply(ind[, , i], 2, sum)
  }
  
  # this gives the total number of frequencies used in all predictors for each central frequency
  totFreqByCol <- apply(ind, 2, sum)
  
  # hFreq <- out$offsets[hIdx]
  
  hOut <- .Fortran("calltf"
                   , d1 = as.double(d1)
                   , d2 = as.double(d2)
                   , ndata = as.integer(ndata)
                   , ndata2 = as.integer(ndata2)
                   , npred = as.integer(dim(d2)[2])
                   , block_size = as.integer(blockSize)
                   , block_size2 = as.integer(blockSize2)
                   , overlap = as.double(overlap)
                   , dt = as.double(dt)
                   , dt2 = as.double(dt2)
                   , nw = as.double(nw)
                   , nw2 = as.double(nw2)
                   , k = as.integer(k)
                   , nFFT = as.integer(nFFT)
                   , nFFT2 = as.integer(nFFT2)
                   , fRatio = as.integer(fRatio)
                   , freq_range_idx = as.integer(freqRangeIdx)
                   , max_freq_offset_idx = as.integer(maxOffIdx)
                   , H = complex(numCol * totalOffsets)
                   , coh_nrow = as.integer(numRow)
                   , coh_ncol = as.integer(numCol)
                   , totFreqByCol = as.integer(totFreqByCol)
                   , totFreqByColByPred = as.integer(totFreqByColByPred)
                   , total_offsets = as.integer(length(hIdx))
                   , col_order = as.integer(ind)
                   , hPredBreak = as.integer(hPredBreak)
                   , hIdx = as.integer(hIdx)
                   , nhIdx = as.integer(length(hIdx)))
  
  H <- matrix(out$H, nrow = numCol, ncol = length(hIdx))
  
  H
}