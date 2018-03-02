#' Calculate the multitaper coherence estimate
#' 
#' Calculates the magnitude squared coherence estimate
#' 
#' @param d1 A \code{numeric} vector containing the first series (central frequency is used in this series).
#' @param d2 A \code{numeric} vector containing the second series (same length as \code{d1}).
#' @param ndata The total length of d1.
#' @param ndata The total length of d2.
#' @param blockSize The length of a single block in series d1 to use (if blocking).
#' @param blockSize2 The length of a single block in series d2 to use (if blocking).
#' @param overlap A \code{numeric} value in the range [0, 1) indicating the proporation of 
#' overlap between neighbouring blocks.
#' @param dt The sampling of rate of series d1 in seconds.
#' @param dt2 The sampling of rate of series d2 in seconds.
#' @param nw time-bandwidth parameter for multitaper (series d1)
#' @param nw2 time-bandwidth parameter for multitaper (series d2)
#' @param k number of tapers to use (k < 2*nw)
#' @param nFFT the number of frequency bins to use for series d1 (nFFT > 2*ndata)
#' @param nFFT2 the number of frequency bins to use for series d2 (nFFT2 > 2*ndata2)
#' @param freqRange A vector with 2 elements containing the start and end frequencies (in Hz) 
#' over which to calculate the coherence.
#' @param maxFreqOffset Every pair of frequencies between f1 (series 1) 
#' and f1 +/- maxFreqOffset (series 2) will be calculated (f1 + maxFreqOffset < nyquist)
#' @param calcType An \code{integer} value indicating how averaging over blocks should 
#' be performed:
#' 1 - calculate the MSC on each block, then average;
#' 2 - calculate the cross and auto spectra on each block, average each quantity 
#' across blocks, then calculate the coherency;
#' 3 - calculate the coherency on each block, then average
#' 4 - returns the minimum MSC across blocks
#' @param forward An \code{integer} indicating whether the forward (1) or reverse (0) coherence
#' should be calculated.
#' @param conv_msc2norm An \code{integer} indicating whether the MSC should be converted 
#' to a Normal distribution (1) or returned as the MSC (0).
#' @param name1 a \code{character} string giving the name of the first data set used.
#' @param name2 a \code{character} string giving the name of the second data set used.
#' 
#' @details Things to think about: 1) blockSize and overlap are going to be very important considerations... 
#' Need to be careful that "block_incr" in mtm_mod.f95 is (are) being calculated correctly.
#' I.e., YOU, the user need to be cogniscent of the values you're using -_- .
#' 
#' If you don't set nw2, this function makes the effective w's the same (i.e., B1 = B2 in BT = NW)
#' 
#' @export
#' @useDynLib transfer2
coherence <- function(d1, d2, ndata = length(d1), ndata2 = length(d2)
                      , blockSize = ndata, blockSize2 = ndata2, overlap = 0
                      , dt = 1, dt2 = dt, nw = 4, nw2 = NULL, k = 7
                      , nFFT = NULL, nFFT2 = nFFT
                      , freqRange = NULL
                      , maxFreqOffset = 0, calcType = 1, forward = 1
                      , conv_msc2norm = 0
                      , name1 = "d1", name2 = "d2")
{
  # You need to be hella careful with nFFT and nFFT2:
  # nFFT2 / nFFT needs to be an integer, otherwise... things will go poorly
  # One (me) assumes that you will be using powers of 2 ... 
  if (dt < dt2){
    stop("d1 is the central frequency series.  You should make the faster sampled series d2.")
  }
  
  if (is.null(nFFT) || nFFT < blockSize) {
    nFFT <- 2^(floor(log2(blockSize))+2)
  }
  
  dtRatio <- dt1 / dt2
  
  if (is.null(nFFT2) || nFFT2 < blockSize2){
    nFFT2 <- nFFT * dtRatio #2^(floor(log2(blockSize2))+2)
  }
  
  fRatio <- (nFFT2 * dt2) / (nFFT * dt)
  
  if ( (fRatio - trunc(fRatio) != 0) | (dtRatio - trunc(dtRatio) != 0) ){
    stop("You should ensure that your choices of nFFT and nFFT2 work with dt and dt2.")
  }
  
  # if nw2 is null, make the effective W's (B - has units) the same
  if (is.null(nw2)){
    # w1 <- nw / (blockSize*dt) 
    # nw2 <- w1*blockSize2*dt2
    # T1*B = NW --> dt1*N1*B = NW1 --> B = NW / (dt1*N1) --> NW2 = dt2*N2*B = [(dt2*N2)/(dt1*N1)] * nw 
    # nw2 <- (dt2*blockSize2*nw) / (dt*blockSize)
    determineNW2(n1 = blockSize, n2 = blockSize2, dt1 = dt, dt2 = dt2, nw1 = nw)
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
  
  nrow = 2*maxOffIdx + 1
  ncol = freqRangeIdx[2] - freqRangeIdx[1] + 1
  
  # print(maxFreqOffset)
  
  # print(paste0("# offsets: ", nrow, " and # freqs: ", ncol))
  # .Fortran("dpss", as.integer(ndata), as.integer(k), as.double(nw), double(ndata*k), double(k))  
  out <- .Fortran("callblockcoh", d1 = as.double(d1), d2 = as.double(d2)
                  , ndata = as.integer(ndata), ndata2 = as.integer(ndata2)
                  , block_size = as.integer(blockSize), block_size2 = as.integer(blockSize2)
                  , overlap = as.double(overlap), dt = as.double(dt), dt2 = as.double(dt2)
                  , nw = as.double(nw), nw2 = as.double(nw2), k = as.integer(k)
                  , nFFT = as.integer(nFFT), nFFT2 = as.integer(nFFT2)
                  , fRatio = as.integer(fRatio)
                  , coh = complex(nrow*ncol), cohnrow = as.integer(nrow)
                  , cohncol = as.integer(ncol)
                  , freq = as.double(freq)
                  , offsets = double(nrow), freq_range_idx = as.integer(freqRangeIdx)
                  , max_freq_offset_idx = as.integer(maxOffIdx)
                  , calc_type = as.integer(calcType)
                  , conv_msc2norm = as.integer(conv_msc2norm)
                  , is_forward = as.integer(forward))
  
  if (calcType == 1){
    calcTypeDesc <- "Average of MSCs"
  } else if (calcType == 2){
    calcTypeDesc <- "MSC of Averages"
  } else if (calcType == 3){
    calcTypeDesc <- "MSC of Average Coherency"
  } else if (calcType == 4){
    calcTypeDesc <- "Minimum of MSCs"
  }
  
  # provides all the argument info used to calculate the coherence
  info <- list(d1 = name1, d2 = name2, ndata = ndata, ndata2 = ndata2
               , blockSize = blockSize, blockSize2 = blockSize2
               , overlap = overlap, nblocks = ndata / (blockSize * (1 - overlap))
               , deltat = dt, deltat2 = dt2
               , nw = nw, k = k, nFFT = nFFT, freqRange = freqRange
               , freqRangeIdx = freqRangeIdx, maxFreqOffset = maxFreqOffset
               , maxFreqOffsetIdx = maxOffIdx
               , calcType = calcType, calcTypeDesc = calcTypeDesc, forward = forward
               , conv_msc2norm = conv_msc2norm)
  
  # list(coh = matrix(out$coh, nrow = nrow))
  if (calcType == 1 | calcType == 4){
    list(coh = matrix(Re(out$coh), nrow = nrow, ncol = ncol), offset = out$offsets
         , bandfreq = freq[freqRangeIdx[1]:freqRangeIdx[2]]
         , info = info)
  } else {
    list(coh = matrix(out$coh, nrow = nrow, ncol = ncol), offset = out$offsets
         , bandfreq = freq[freqRangeIdx[1]:freqRangeIdx[2]]
         , info = info)
  }
}


#' Transform coherence to Guassian distribution
#' 
#' Converts estimates of MSC to an almost Gaussian distribution
#' 
#' @param c a real-valued vector of the coherence
#' @param k the number of tapers used
#' @param msc a \code{logical} indicating whether the coherence passed was actually the MSC
#' 
#' @export
mscQTransform <- function(c, k, msc = TRUE){
  if (msc){
    sqrt(2*k-2) * atanh(sqrt(c))
  } else {
    sqrt(2*k-2) * atanh(c)
  }
}

#' Normal Transform for Magnitude Squared Coherence
#' 
#' Converts the MSC to a Gaussian distribution
#' 
#' @param msc A vector containing the magnitude squared coherence.
#' @param dof Degrees of freedom of the msc
#' 
#' @details This code was provided by Emily Somerset.
#' 
#' @export
msc2norm <- function(msc, dof){
  tran <- 1 - (1-msc)^(dof-1) #
  erf.inv <- function(x) { qnorm((x + 1)/2)/sqrt(2) }
  
  0 + sqrt(2)*1*erf.inv(2*tran - 1)
}

#' Determines the local maxes of MSC
#' 
#' Finds the local max of the normal transformed MSC by checking the neighbours.
#' 
#' @param msc the normal transformed magnitude square coherence (see \link{mscQTransform})
#' @param k an \code{integer} indicating the number of tapers used to calculate the MSC
#' @param cutoff a value in [0, 1] representing the significance level.
#' 
#' @export
findLocalMscMax <- function(msc, k, cutoff){
  lev <- 1 - (1-cutoff)^(1/(k-1))
  maxInd <- which(msc > lev)## based on normal transform.. qnorm(cutoff, mean = 1/sqrt(2*k-2), 1))
  maxes <- c()
  
  if (length(maxInd) == 0){
    return(maxes)
  }
  
  for (i in 1:length(maxInd)){
    if (maxInd[i] == 1 || maxInd[i] == length(msc)){
      next
    }
    
    if (msc[maxInd[i]] > msc[maxInd[i]-1] && 
        msc[maxInd[i]] > msc[maxInd[i]+1]){
      maxes <- c(maxes, maxInd[i])
    }
  }
  
  maxes
}

#' Determines the local maxes of F-test values
#' 
#' Finds the local max of the multitaper harmonic F-test by checking the neighbours.
#' 
#' @param fval a \code{vector} containing the value sof the f-test.
#' @param k an \code{integer} indicating the number of tapers used to calculate the F-values.
#' @param cutoff a value in [0, 1] representing the significance level.
#' 
#' @export
findLocalFvalMax <- function(fval, k, cutoff){
  maxInd <- which(fval > qf(cutoff, 2, 2*k-2))
  maxes <- c()
  
  if (length(maxInd) == 0){
    return(maxes)
  }
  
  for (i in 1:length(maxInd)){
    if (maxInd[i] == 1 || maxInd[i] == length(fval)){
      next
    }
    
    if (fval[maxInd[i]] > fval[maxInd[i]-1] && 
        fval[maxInd[i]] > fval[maxInd[i]+1]){
      maxes <- c(maxes, maxInd[i])
    }
  }
  
  maxes
}
