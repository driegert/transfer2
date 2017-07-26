#' @export
#' @useDynLib transfer2
fortTest <- function(x){
  message("Holy poop-farts!")
  out <- .Fortran("adder", as.double(x))
  out
}

#' @export
#' @useDynLib transfer2
forTest2 <- function(x){
  out <- .Fortran("test2", as.integer(x))
  out
}

#' @export
fortTest3 <- function(){
  out <- .Fortran("test3", as.integer(3))
  out
}

#' @export
#' @useDynLib transfer2
testmodule <- function(x){
  out <- .Fortran("testem", x = as.integer(x))
  out
}

#' @export
# #' @useDynLib libdlr
testDPSS <- function(n, nw, k){
  out <- .Fortran("dpss", n = as.integer(n), k = as.integer(k), nw = as.double(nw)
                  , v = double(n*k), eigen = double(k))
  out
}

#' @export
testcd <- function(){
  cs <- matrix(complex(real = 1, imaginary = 1), ncol = 10, nrow = 11)
  s1 <- 1:129; s2 = 1:129

  out <- .Fortran("testcd", cs = as.complex(cs), s1 = as.double(s1), s2 = as.double(s2)
                  , cnrow = as.integer(11), cncol = as.integer(10)
                  , idx_start = as.integer(51), idx_end = as.integer(60)
                  , max_idx = as.integer(5), nFFT = as.integer(256))
  out
}

#' @export
testdpss <- function(n, nw, k){
  nFFT <- 2^(floor(log2(n)) + 2)
  out <- .Fortran("testdpss", n = as.integer(n), nw = as.double(nw)
                  , k = as.integer(k), vec = double(n*k)
                  , val = double(k), nFFT = as.integer(32))
  out
}

#' @export
testratios <- function(){
  ndata = 16; nw = 4; k = 7; nFFT = 32
  out <- .Fortran("testratios", ndata = as.integer(ndata), nw = as.double(nw), k = as.integer(k)
                  , ratios = double(ndata)
                  , vec = double(ndata*k), val = double(k), nFFT = as.integer(nFFT)
                  , cx = complex(nFFT*k))
  list(ndata = out$ndata, nw = out$nw, k = out$k, ratios = out$ratios
       , v = matrix(out$vec, nrow = ndata, ncol = k), ev = out$val
       , nFFT = out$nFFT, cx = matrix(Re(out$cx), nrow = nFFT, ncol = k))
}

#' @export
testeigen <- function(d, nw, k, nFFT, dt){
  n <- length(d)
  nfreq <- nFFT/2+1
  out <- .Fortran("testeigen", d = as.double(d)
                  , ndata = as.integer(n)
                  , nw = as.double(nw), k = as.integer(k)
                  , yk = complex(nfreq*k)
                  , nFFT = as.integer(nFFT)
                  , dt = as.double(dt))
  list(yk = matrix(out$yk, ncol = k, nrow = nfreq))
}

#' @export
testweights <- function(dat, yk, k, dt, var, nFFT, nw, eval){
  n <- length(dat)
  nfreq <- nFFT/2+1
  var <- var * (n-1)/n
  out <- .Fortran("testweights", d = as.double(dat), yk = as.complex(yk)
                  , dk = double(nfreq * k), k = as.integer(k)
                  , var = as.double(var), dt = as.double(dt)
                  , ndata = as.integer(n)
                  , nFFT = as.integer(nFFT)
                  , nw = as.double(nw)
                  , eval = as.double(eval))
  out
}

#' @export
weightwrapper <- function(yk, nFFT, k, ev, var, dt, n){
  sa = abs(yk)^2
  nfreq = nFFT/2+1
  var <- var * (n-1)/n
  out <- .Fortran("weightwrapper", sa = as.double(sa)
                  , dk = double(nfreq*k), nFFT = as.integer(nFFT)
                  , k = as.integer(k), ev = as.double(ev)
                  , var = as.double(var), dt = as.double(dt))
  
  out
}

#' @export
testfft <- function(dat, nFFT){
  n <- length(dat)
  out <- .Fortran("testfft", d = as.double(dat), ndata = as.integer(n)
                  , nFFT = as.integer(nFFT), c = complex(nFFT))
  out
}

#' @export
testpfft <- function(dat, nFFT){
  n <- length(dat)
  out <- .Fortran("testpfft", d = as.double(dat), ndata = as.integer(n)
                  , nFFT = as.integer(nFFT), c = complex(nFFT))
  out
}

#' @export
testcross <- function(yk1, yk2, nFFT, idx_start, idx_end, idx_max){
  k <- dim(yk1)[2]
  nfreq <- nFFT/2+1
  crossSize <- c(2*idx_max + 1, idx_end - idx_start + 1)
  out <- .Fortran("testcross", yk1 = as.complex(yk1), yk2 = as.complex(yk2)
                  , cs12 = complex(crossSize[1] * crossSize[2]), k = as.integer(k)
                  , nFFT = as.integer(nFFT), idx_start = as.integer(idx_start)
                  , idx_end = as.integer(idx_end), idx_max = as.integer(idx_max))
  out
}

#' @export
testcross2 <- function(d1, d2, nw, k, nFFT, idx_start, idx_end, idx_max, dt){
  n <- length(d1)
  nfreq <- nFFT/2+1
  crossSize <- c(2*idx_max + 1, idx_end - idx_start + 1)
  out <- .Fortran("testcross2", d1 = as.double(d1), d2 = as.double(d2)
                  , ndata = as.integer(n), cs12 = complex(crossSize[1] * crossSize[2])
                  , nw = as.double(nw), k = as.integer(k), nFFT = as.integer(nFFT)
                  , idx_start = as.integer(idx_start), idx_end = as.integer(idx_end)
                  , idx_max = as.integer(idx_max), dt = as.double(dt))
  out
}

#' Test out to see if the regression stuff works... 
#' @export
testSvdReg <- function(Y, X){
  m <- dim(X)[1]
  n <- dim(X)[2]
  
  out <- .Fortran("svdRegTest", Y = as.complex(Y), X = as.complex(X), m = as.integer(m)
                  , n = as.integer(n), beta = complex(n), stdErr = double(n), svd_ev = double(n))
  out
}

#' @export
testZsvd <- function(Y, X){
  m <- dim(X)[1]
  n <- dim(X)[2]
  
  out <- .Fortran("zsvd", Y = as.complex(Y), X = as.complex(X), m = as.integer(m), n = as.integer(n)
                  ,  u = complex(m*n), s = double(n), vt = complex(n*n))
  out
}