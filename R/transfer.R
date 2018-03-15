#' Estimates the transfer functions
#' 
#' Estimates a transfer function using frequency offsets between a response and multiple inputs
#' 
#' @param d1 A single column \code{data.frame} containing the response.
#' @param d2 A \code{data.frame} whose columns contain the predictors.
#' @param blockSize The number of points per block to use for the response (default = dim(d1)[1])
#' @param blockSize2 The number of points per block to use for the inputs (default = dim(d2)[2])
#' @param overlap A number between 0 and 1 representing the proportion of overlap between adjacent
#' blocks.
#' @param dt A value representing the sampling rate of the response (in seconds).
#' @param dt2 A value representing the sampling rate of the predictors (in seconds).
#' @param nw The time bandwidth parameter to use with the response.
#' @param nw2 The time bandwidth parameter to use with the predictors (only set this if you 
#' know what you're doing -_- )
#' @param k The number of tapers to use (the same for response and predictors).
#' @param nFFT The number of frequency bins to be used for the response (sets zero-padding).
#' @param nFFT2 The number of frequency bins to be used for predictors (sets zero-padding) - 
#' (only set this if you know what you're doing).
#' @param freqRange A vector of 2 values containing the start and end values of the frequency band 
#' for the response, over which to estimate the transfer function (in Hz).
#' @param maxFreqOffset The maximum from the central frequency to use when estimating the transfer 
#' functions (careful here too).
#' @param nOff An \code{integer} representing the number of offset frequencies that are allowed to be used.
#' @param forceZeroFreq A \code{logical} indicating whether to always use the zero-offset frequency
#' (default = TRUE).
#' @param sigLevel A value between 0 and 1 representing the confidence level to use when determining whether 
#' the coherence between response and predictor should be considered statistically significant.
#' @param standardize A \code{logical} indicating whether the data should have its mean subtracted off and 
#' standard deviation divided out.
#' @param name1 I don't think this is used ... 
#' @param name2 I don't think this is used ... yet.
#' 
#' @export
#' @useDynLib transfer2
tf <- function(d1, d2
               , blockSize = dim(d1)[1], blockSize2 = dim(d2)[1], overlap = 0
               , dt = 1, dt2 = dt, nw = 4, nw2 = NULL, k = 7
               , nFFT = NULL, nFFT2 = NULL
               , freqRange = NULL, maxFreqOffset = 0, nOff = -1
               , forceZeroFreq = TRUE
               , sigLevel = 0.99
               , standardize = TRUE
               , name1 = names(d1), name2 = names(d2)){
  
  ndata <- dim(d1)[1]
  ndata2 <- dim(d2)[1]
  npred <- dim(d2)[2]
  
  # You need to be hella careful with nFFT and nFFT2:
  # nFFT2 / nFFT needs to be an integer, otherwise... things will go poorly
  # One (me) assumes that you will be using powers of 2 ... 
  if (dt < dt2){
    stop("d1 is the central frequency series.  You should make the faster sampled series d2.")
  }
  
  if (standardize){
    d1 <- as.data.frame(lapply(d1, function(x){ (x - mean(x)) / sd(x) }))
    d2 <- as.data.frame(lapply(d2, function(x){ (x - mean(x)) / sd(x) }))
  }
  
  if (is.null(nFFT) || nFFT < blockSize) {
    nFFT <- 2^(floor(log2(blockSize))+3)
  }
  
  dtRatio <- dt / dt2
  
  if (is.null(nFFT2) || nFFT2 < blockSize2){
    nFFT2 <- nFFT * dtRatio # 2^(floor(log2(blockSize2))+3)
  }
  
  # M1 / (1/dt1) == a * M2 / (1/dt2), solve for 'a'
  fRatio <- determineFreqRatio(dt, dt2, nFFT, nFFT2)#(nFFT2 * dt2) / (nFFT * dt)
  
  if (fRatio - trunc(fRatio) != 0){
    stop("You should ensure that your choices of nFFT and nFFT2 work with dt and dt2.\n 
         i.e., df2 / df1 should be an integer.")
  }
  
  # if nw2 is null, make the effective W's the same
  w1 <- nw / (dt*blockSize) # effective W -> needed for the deadzone
  if (is.null(nw2)){
    # nw2 <- w1*blockSize2
    nw2 <- determineNW2(n1 = blockSize, n2 = blockSize2, dt1 = dt, dt2 = dt2, nw1 = nw)
  }
  
  df <- 1 / (dt*nFFT)
  
  # this if-statement is here because there were issues when freqRange == NULL
  if (is.null(freqRange)){
    warnings("freqRange == NULL: Setting freqRange to full positive band.")
    freqRange <- c(0, 1/(2*dt))
    # maxFreqOffset <- 0
    
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
  
  # number of predictors has to be <= number of samples.
  if ((forceZeroFreq && (npred*(nOff-1) >= (nblocks*k))) || (forceZeroFreq && (npred*nOff >= (nblocks*k)))){
    stop("Set nOffs so that the number of 'samples', i.e., blocks * tapers is less than the potential number 
         of predictors, i.e., nOff + zeroFreq (if forceZeroFreq == TRUE).\n
         Number of coefficients estimated must be less than (or equal, but I don't allow that -_-) to the number 
         of samples available.")
  }
  
  # determine which frequencies get used
  level <- qnorm(sigLevel) / sqrt(nblocks)
  
  # Just need to do a "normal" transfer function regression - write some code in Fortran.
  if (maxOffIdx == 0){
    out <- .Fortran("calltfzero"
                    , d1 = as.double(d1[,1])
                    , d2 = as.double(as.matrix(d2))
                    , ndata = as.integer(ndata)
                    , ndata2 = as.integer(ndata2)
                    , npred = as.integer(npred)
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
                    , H = complex(numCol * npred)
                    , n_row_H = as.integer(numCol))
    
    H <- matrix(out$H, nrow = numCol, ncol = npred)
    hIdx <- rep(1, npred)
    hPredBreak <- 1:4
  } else {
    if (forceZeroFreq & (maxOffIdx < zeroDeadZone)){
      stop("maxFreqOffset should be larger than (B1)/3 in order to use the zero offset frequency.")
    }
    
    # calculate the pairwise MSC's
    ind <- array(data = 0, dim = c(numRow, numCol, ncol(d2)))
    for (i in 1:ncol(d2)){
      out <- .Fortran("callblockcoh"
                      , d1 = as.double(d1[, 1])
                      , d2 = as.double(d2[, i])
                      , ndata = as.integer(ndata)
                      , ndata2 = as.integer(ndata2)
                      , block_size = as.integer(blockSize)
                      , block_size2 = as.integer(blockSize2)
                      , overlap = as.double(overlap)
                      , dt = as.double(dt)
                      , dt2 = as.double(dt2)
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
    totalOffsets <- length(hIdx)
    # browser()
    hOut <- .Fortran("calltf"
                     , d1 = as.double(d1[, 1])
                     , d2 = as.double(as.matrix(d2))
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
                     , total_offsets = as.integer(totalOffsets)
                     , col_order = as.integer(ind)
                     , hPredBreak = as.integer(hPredBreak)
                     , hIdx = as.integer(hIdx)
                     , nhIdx = as.integer(length(hIdx)))
    H <- matrix(hOut$H, nrow = numCol, ncol = totalOffsets)
  }
  
  # we need to obtain the value of the offsets to be used - probably actualy index and also in terms of frequency
  offIdxFromCent <- seq(-(fRatio*maxOffIdx), fRatio*maxOffIdx, by = fRatio)
  offFreqFromCent <- offIdxFromCent * df
  
  # create a data.frame of info for the transfer function matrix.
  for (i in 1:ncol(d2)){
    if (i == 1){
      Hinfo <- data.frame(predictor = as.character(names(d2)[i])
                          , Hcolumn = hPredBreak[i]:(hPredBreak[i+1]-1)
                          , freqOffset = offFreqFromCent[hIdx[ hPredBreak[i]:(hPredBreak[i+1] - 1) ]]
                          , idxOffset = offIdxFromCent[hIdx[ hPredBreak[i]:(hPredBreak[i+1] - 1) ]])
    } else {
      Hinfo <- rbind(Hinfo,data.frame(predictor = as.character(names(d2)[i])
                                      , Hcolumn = hPredBreak[i]:(hPredBreak[i+1]-1)
                                      , freqOffset = offFreqFromCent[hIdx[ hPredBreak[i]:(hPredBreak[i+1] - 1) ]]
                                      , idxOffset = offIdxFromCent[hIdx[ hPredBreak[i]:(hPredBreak[i+1] - 1) ]])
                     )
    }
    # hColNames <- c(hColNames
    #                , paste0(names(d2)[i], "..", offFreqFromCent[hIdx[ hPredBreak[i]:(hPredBreak[i+1] - 1) ]]))
  }
  
  info <- list(ndata = ndata, ndata2 = ndata2, blockSize = blockSize, blockSize2 = blockSize2
               , dt = dt, dt2 = dt2, nw = nw, nw2 = nw2, k = k, nFFT = nFFT, nFFT2 = nFFT2
               , freqRange = freqRange, freqRangeIdx = freqRangeIdx
               , maxFreqOffset = maxFreqOffset, maxOffsetIdx = maxOffIdx
               , nOff = nOff, forceZeroFreq = forceZeroFreq, sigLevel = sigLevel
               , standardize = standardize)
  
  list(H = H, Hinfo = Hinfo, info = info)
}

### maybe want this later?  not sure if the eigencoefficients are quite what we want.
ykPredict <- function(H, Hinfo, info, d2){
  # make sure column names of d2 match those in Hinfo
  ###
  ###
  spec <- list()
  yk2 <- list()
  for (i in 1:ncol(d2)){
    spec[[i]] <- multitaper::spec.mtm(d2[, i], nw = info$nw2, k = info$k, deltat = info$dt2
                                      , nFFT = info$nFFT2
                                      , returnInternals = TRUE, plot = FALSE)
    yk2[[i]] <- spec[[i]]$mtm$eigenCoefs * spec[[i]]$mtm$eigenCoefWt
  }
}

#' Spectrum prediction
#' Predicts the spectrum based on an estimated transfer function and new data provided.
#' 
#' @param H The H returned by tf().
#' @param Hinfo The Hinfo returned by tf().
#' @param info The info returned by tf().
#' @param d2 A \code{data.frame} containing the new data.  Must have the same column names as the original 
#' data used in the transfer function estimation.
#' 
#' @export
specPredict <- function(H, Hinfo, info, d2){
  predNames <- names(d2)
  spec <- list()
  for (i in 1:length(predNames)){
    spec[[predNames[i]]] <- multitaper::spec.mtm(d2[, i], nw = info$nw2, k = info$k, deltat = info$dt2
                                                 , nFFT = info$nFFT2
                                                 , center = 'none'
                                                 , returnInternals = TRUE, plot = FALSE)$spec
  }
  
  fullFreqRange <- info$freqRangeIdx[1]:info$freqRangeIdx[2]
  
  sRecon <- rep(0, length(fullFreqRange))
  
  for (i in 1:dim(Hinfo)[1]){
    offIdx <- fullFreqRange + Hinfo$idxOffset[i]
    offIdx[offIdx <= 0] <- abs(offIdx[offIdx <= 0]) + 2
    sRecon <- sRecon + (abs(H[, i])^2) * spec[[ Hinfo$predictor[i] ]][offIdx]
  }
  
  sRecon
}