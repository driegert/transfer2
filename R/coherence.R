#' Calculate the multitaper coherence estimate
#' 
#' Calculates the magnitude squared coherence estimate
#' 
#' @param d1 A \code{numeric} vector containing the first series.
#' @param d2 A \code{numeric} vector containing the second series (same length as \code{d1}).
#' @param ndata The total length of d1 or d2.
#' @param blockSize The length of a single block to use (if blocking)
#' @param overlap A \code{numeric} value in the range [0, 1) indicating the proporation of 
#' overlap between neighbouring blocks.
#' @param dt The sampling rate in seconds.
#' @param nw time-bandwidth parameter for multitaper
#' @param k number of tapers to use (k < 2*nw)
#' @param nFFT the number of frequency bins to use (nFFT > 2*ndata)
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
#' @param forward An \code{integer} indicating whether the forward (1) or reverse (0) coherence
#' should be calculated.
#' 
#' @export
#' @useDynLib transfer2
coherence <- function(d1, d2, ndata = length(d1), blockSize = ndata, overlap = 0
                      , dt = 1, nw = 4, k = 7, nFFT = NULL, freqRange = NULL
                      , maxFreqOffset = 0, calcType = 1, forward = 1)
{
  if (is.null(nFFT) || nFFT < blockSize) {
    nFFT <- 2^(floor(log2(blockSize))+2)
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
                  , ndata = as.integer(ndata), block_size = as.integer(blockSize)
                  , overlap = as.double(overlap), dt = as.double(dt)
                  , nw = as.double(nw), k = as.integer(k), nFFT = as.integer(nFFT)
                  , coh = complex(nrow*ncol), cohnrow = as.integer(nrow)
                  , cohncol = as.integer(ncol)
                  , freq = as.double(freq)
                  , offsets = double(nrow), freq_range_idx = as.integer(freqRangeIdx)
                  , max_freq_offset_idx = as.integer(maxOffIdx)
                  , calc_type = as.integer(calcType)
                  , is_forward = as.integer(forward))
  
  # list(coh = matrix(out$coh, nrow = nrow))
  if (calcType == 1){
    list(coh = matrix(Re(out$coh), nrow = nrow, ncol = ncol), offset = out$offsets
         , bandfreq = freq[freqRangeIdx[1]:freqRangeIdx[2]])
  } else {
    list(coh = matrix(out$coh, nrow = nrow, ncol = ncol), offset = out$offsets
         , bandfreq = freq[freqRangeIdx[1]:freqRangeIdx[2]])
  }
}

#' Transform coherence to Guassian distribution
#' 
#' Converts estimates of MSC to an almost Gaussian distribution
#' 
#' @param c a real-valued vector of the coherence
#' @param k the number of tapers used
#' @parma msc a \code{logical} indicating whether the coherence passed was actually the MSC
#' 
#' @export
mscQTransform <- function(c, k, msc = TRUE){
  if (msc){
    sqrt(2*k-2) * atanh(sqrt(c))
  } else {
    sqrt(2*k-2) * atanh(c)
  }
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
