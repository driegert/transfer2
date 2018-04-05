# Determines nw2 such that the effective w1 and w2 are equal.
# This should mean that nw2 is at least as large as nw1
determineNW2 <- function(n1, n2, dt1, dt2, nw1){
  (dt2*n2*nw1) / (dt1*n1)
}

# Provides the ratio 
determineFreqRatio <- function(dt1, dt2, nFFT1, nFFT2){
  (nFFT2 * dt2) / (nFFT1 * dt1)
}

#' @export
df.std <- function(x){
  stopifnot(class(x) == "data.frame")
  
  as.data.frame(lapply(x, function(z){ (z - mean(z)) / sd(z) }))
}

#' Determines Offset Frequencies Used
#' TESTING - NOT FINISHED.
offsetFrequencies <- function(){
  df <- 1 / (dt*nFFT)
  
  maxOffIdx <- floor(maxFreqOffset / df)
}


#' Convolve Complex-Valued Impulse Response
#' 
#' @param x A \code{vector} on which to apply the impulse response.
#' @param filter A \code{vector} the impulse response, must be odd-length (for now).
#' 
#' @details Taken from transfer package (Aaron and I wrote that).
#' This only accepts two-sided filtering for now.
#' 
#' @export
zFilter <- function (x, filter, realPart = TRUE) 
{
  N <- length(x)
  nFilt <- length(filter)
  if (nFilt > N) {
    stop("Filter is longer than the time series.")
  }
  if (nFilt%%2 == 0) {
    stop("Filter length should be odd ... or modify this function to work with even length.")
  }
  M <- 2^(floor(log2(N)) + 3)
  result <- fft(fft(c(x, rep(0, M - N))) * fft(c(filter, rep(0, 
                                                             M - nFilt))), inverse = TRUE)[1:N]/M
  result[1:(nFilt - 1)] <- NA
  if (realPart){
    Re(result)
  } else {
    result
  }
}

#' I don't know if I quite trust this.. worth another look
#' @export
H2zero <- function(H, freqBand){
  freqBandIdx <- c(max(1, floor(freqBand[1] / H$info$df)), ceiling(freqBand[2] / H$info$df))
  
  if (any((H$info$freqRangeIdx - freqBandIdx) < 0)){
    stop("freqBand must be a subset of H$info$freqRange")
  }
  
  freq <- H$info$freqRangeIdx[1]:H$info$freqRangeIdx[2]
  
  idx <- which(freq == freqBandIdx[1]):which(freq == freqBandIdx[2])
  
  H$H[idx, ] <- complex(real = 0, imaginary = 0)
  
  H
}

#' this only works for a very specific situation
#' @export
fixDeadBand <- function(msc, zeroOffsetIdx, freqRangeIdx, band, df){
  # if ((band[1] > 0 | band[2] < 0) | (freqRangeIdx[1] != 1)){
  #   stop("fixDeadBand doesn't work unless you start at 0 and you're taking a band around 0.")
  # }
  
  ## assuming symmetry around f = 0:
  nIdx <- ceiling(band / df)
  
  top <- zeroOffsetIdx + nIdx
  bottom <- zeroOffsetIdx - nIdx
  
  # everything is backwards - we need the reverse diagonals.  Also, in fields::image.plot(), the
  # y-axis is "reversed" as it were... so thinking about this based on those plots is difficult..
  matTop <- lower.tri(matrix(0, top, top), diag = TRUE)
  matBot <- upper.tri(matrix(0, nrow = bottom, ncol = top), diag = FALSE)
  matTop[(top-bottom+1):top, 1:top] <- matBot & matTop[(top-bottom+1):top, 1:top]
  
  msc[top:1, 1:top][matTop] <- 0
  
  msc
}