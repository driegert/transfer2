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
#' @param deadBand A 2-vector containing the start and end region of a band to ignore when estimating the 
#' transfer functions.  This would occur due to high pass filtering for example (low frequency trend removal).
#' 
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
               , deadBand = NULL
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
    d1.stdPar <- data.frame(variable = names(d1), mean = mean(d1[, 1]), sd = sd(d1[, 1]), stringsAsFactors = FALSE)
    d1 <- as.data.frame(lapply(d1, function(x){ (x - mean(x)) / sd(x) }))
    pnam <- names(d2)
    d2.stdPar <- data.frame(variable = pnam, mean = rep(-999.9, ncol(d2)), sd = rep(-999.9, ncol(d2))
                            , stringsAsFactors = FALSE)
    for (i in 1:ncol(d2)){
      d2.stdPar[which(d2.stdPar$variable == pnam[i]), c("mean", "sd")] <- c(mean(d2[, pnam[i]]), sd(d2[, pnam[i]]))
    }
    d2 <- as.data.frame(lapply(d2, function(x){ (x - mean(x)) / sd(x) }))
  } else {
    d1.stdPar <- NULL
    d2.stdPar <- NULL
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
  df2 <- 1 / (dt2 * nFFT2)
  
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
                      , conv_msc2norm = as.integer(0)
                      , is_forward = as.integer(1))

      #########################################################
      # FIX THIS ... error function estimation is BAD in msc2norm() (and in the Fortran code too... )
      # this breaks when the MSC is >= 0.94
      out$coh[which(Re(out$coh) > 0.93)] <- complex(real = 0.93, imaginary = 0) ## this is a fucking awful hack.
      out$coh <- msc2norm(Re(out$coh), dof = k)
      #################################################
      
      ##########
      ### there are issues if you try to run this code after filtering (division by close to 0's)
      # the problem occurs when the the offset frequency used close to the filter "hole" in the spectrum
      # so we modify the coherence by setting it to zero in the bands where we "know" issues will occur.
      ###### Have to come back and fix this...
      if (!is.null(deadBand)){
        out$coh <- fixDeadBand(msc = matrix(Re(out$coh), nrow = numRow, ncol = numCol)
                               , zeroOffsetIdx = maxOffIdx + 1
                               , freqRangeIdx = freqRangeIdx, band = deadBand[2], df = df)
      }
      
      out2 <- .Fortran("callMscIndicator", msc = as.double(Re(out$coh)), nrow = as.integer(numRow)
                       , ncol = as.integer(numCol), ind = integer(numRow * numCol)
                       , level = as.double(level), nOff = as.integer(nOff))
      
      ind[, , i] <- matrix(out2$ind, nrow = numRow, ncol = numCol)
      
      
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
               , df = df, df2 = df2
               , nOff = nOff, forceZeroFreq = forceZeroFreq, sigLevel = sigLevel
               , standardize = standardize, d1.stdPar = d1.stdPar
               , d2.stdPar = d2.stdPar)
  
  list(H = H, Hinfo = Hinfo, info = info)
}

### maybe want this later?  not sure if the eigencoefficients are quite what we want.
# ykPredict <- function(H, Hinfo, info, d2){
#   # make sure column names of d2 match those in Hinfo
#   ###
#   ###
#   spec <- list()
#   yk2 <- list()
#   for (i in 1:ncol(d2)){
#     spec[[i]] <- multitaper::spec.mtm(d2[, i], nw = info$nw2, k = info$k, deltat = info$dt2
#                                       , nFFT = info$nFFT2
#                                       , returnInternals = TRUE, plot = FALSE)
#     yk2[[i]] <- spec[[i]]$mtm$eigenCoefs * spec[[i]]$mtm$eigenCoefWt
#   }
# }

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

#' Impulse response from transfer functions
#' Inverse Fourier transforms the transfer function to give the impulse response
#' 
#' @param H A \code{matrix} where each column contains a transfer function - tf()$H would be appropriate.
#' @param n An \code{integer} indicating the half-length of the impulse response to return, L = 2*n+1. 
#' Normally this would be the \code{blockSize2} agrument used in tf().
#' @param realPart A \code{logical} indicating whether to return only the real part of the impulse response
#' (default = TRUE) or return the complex-valued impulse response (FALSE)
#' 
#' @export
# Careful here... we're using R's FFT on the H's calculated in Fortran with FFTpack (or whatever it's called)
ir <- function(H, n, realPart = TRUE){
  stopifnot(is.matrix((H)) | n > ncol(H))
  # "negative" frequencies in the top half of the array
  # ^^ these are conjugated first and reversed (should be conjugate symmetric after)
  Hfull <- rbind(H, Conj(H[(nrow(H) - 1):2, , drop = FALSE]))
  
  # should be real-valued with real-valued data as inputs.
  if (realPart){
    h <- Re(mvfft(Hfull, inverse = TRUE) / nrow(Hfull))
  } else {
    h <- mvfft(Hfull, inverse = TRUE) / nrow(Hfull)
  }
  
  n <- n+1
  # put the impulse response in the correct place in the array in order to use in a convolution - i.e., filter()
  ind <- c((nrow(h)-n+2):nrow(h), 1:n)
  
  h[ind, , drop = FALSE]
}


#' Predict the Time Series from the predictors and impulse response
#' 
#' @param newdata A \code{data.frame} containing the predictors.  Column names must match what is 
#' present in Hinfo.
#' @param ir A \code{matrix} containing the impulse responses corresponding to the rows in Hinfo.
#' @param Hinfo A piece of what is returned by tf().
#' @param info Another piece of what is returned by tf().
#' 
#' @export
predictTs <- function(newdata, ir, Hinfo, info){
  N <- nrow(newdata)
  yhat <- rep(0, N)
  nFlt <- nrow(ir)
  nTrim <- (nFlt - 1) / 2
  
  # standardize the input if necessary.
  if (info$standardize){
    newdata <- df.std(newdata)
  }
  
  for (i in 1:nrow(Hinfo)){
    ## Should this be 2*pi*(1:N) or 2*pi*(0:N) ?
    ## is the negative necessary? I guess it doesn't matter since cos(x) is an even function... 
    tmp <- zFilter(2*cos(-2*pi*(0:(N-1))*(Hinfo$idxOffset[i]/info$nFFT2)) * newdata[, Hinfo$predictor[i]]
                   , ir[, Hinfo$Hcolumn[i]])
    
    yhat <- yhat + tmp[nTrim + 1:N] # this adds NA's at the end, also has NA's at the beginning - I need to look at 
    # zFilter again for how this works by convoling using FFT's ... 
  }
  
  yhat
}