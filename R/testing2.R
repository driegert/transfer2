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