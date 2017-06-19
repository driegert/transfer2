#' @export
#' @useDynLib transfer2
coherence <- function(d1, d2, ndata = length(d1), blockSize = ndata, overlap = 0
                      , dt = 1, nw = 4, k = 7, nFFT = NULL, freqRange = NULL
                      , maxFreqOffset = 0, calcType = 1, forward = 1)
{
  if (is.null(nFFT) || nFFT < blockSize) {
    nFFT <- 2^(floor(log2(blockSize))+2)
  }
  
  if (is.null(freqRange)){
    warnings("freqRange == NULL: Setting freqRange to positive band and maxFreqOffset to 0.")
    freqRange <- c(0, 1/(2*dt))
    maxFreqOffset <- 0
  }
  
  print(paste0("ndata: ", ndata))
  
  df <- 1 / (dt*nFFT)
  maxOffIdx <- ceiling(maxFreqOffset / df)
  freqRangeIdx <- c(floor(freqRange[1] / df), ceiling(freqRange[2]/df))
  freq <- seq(0, 1/(2*dt), 1/(nFFT*dt))
  
  nrow = 2*maxOffIdx + 1
  ncol = freqRangeIdx[2] - freqRangeIdx[1] + 1
  
  # print(maxFreqOffset)
  
  print(paste0("# offsets: ", nrow, " and # freqs: ", ncol))
  
  # .Fortran("dpss", as.integer(ndata), as.integer(k), as.double(nw), double(ndata*k), double(k))  
  out <- .Fortran("callblockcoh", d1 = as.double(d1), d2 = as.double(d2)
                  , ndata = as.integer(ndata), block_size = as.integer(blockSize)
                  , overlap = as.double(overlap), dt = as.double(dt)
                  , nw = as.double(nw), k = as.integer(k), nFFT = as.integer(nFFT)
                  , coh = complex(nrow*ncol), cohnrow = as.integer(nrow)
                  , cohncol = as.integer(ncol)
                  , freq = as.double(freq)
                  , offsets = double(nrow), freq_range = as.double(freqRange)
                  , max_freq_offset = as.double(maxFreqOffset)
                  , calc_type = as.integer(calcType)
                  , is_forward = as.integer(forward))
  
  # list(coh = matrix(out$coh, nrow = nrow))
  list(coh = matrix(out$coh, nrow = nrow, ncol = ncol), offset = out$offsets
       , bandfreq = freq[freqRangeIdx[1]:freqRangeIdx[2]])
}