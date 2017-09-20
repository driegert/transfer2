module mtm_mod
  implicit none
  real*8, parameter :: pi = dble(3.1415926535897932)
  real*8, allocatable :: wsave(:), work(:), v(:, :), ev(:)
  integer :: lenc, lensav, lenwrk, inc, ier, nfreq
  logical :: is_dpss_setup = .false., is_fft_setup = .false.

contains
  subroutine fft_setup(nFFT)
  ! calls the setup routine for fftpack5.1d routines
  ! cfft1f (or cfft1b)
    integer :: nFFT

    lenc = nFFT
    nfreq = nFFT/2 + 1
    inc = 1
    lensav = 2*nFFT + int(log( dble(nFFT) ) / log( dble(2) )) + 4
    lenwrk = 2*nFFT
    allocate(wsave(lensav))
    allocate(work(lenwrk))

    call cfft1i(nFFT, wsave, lensav, ier)
  end subroutine fft_setup

  subroutine fft_cleanup()
  ! Frees up memory again so that this code can be re-called.
    deallocate(wsave, work)
  end subroutine

  subroutine dpss_setup(ndata, nw, k, nFFT)
  ! calculates the dpss tapers and eigenvalues
  !
  ! ndata - integer - length of data to be used
  ! nw    - real*8  - time bandwidth parameter
  ! k     - integer - number of tapers to be used
  ! nFFT  - integer - length to zeropad out to
  !
    integer :: ndata, k, nFFT
    real*8 :: nw
    allocate(v(ndata, k), ev(k))

    call dpss(ndata, k, nw, v, ev)
    call dpssToEigenvalues(ndata, k, nw, nFFT)
  end subroutine dpss_setup

  subroutine dpss_cleanup()
  ! Frees up memory again so that this code can be re-called.
    deallocate(v, ev)
  end subroutine dpss_cleanup

  subroutine cross_spec(yk1, yk2, cs12, k, idx_start, idx_end, max_idx, is_forward)
  ! Calculates the cross spectrum given two eigencoefficient matrices
  !
  ! yk1 - complex*16(nfreq, k) - series 1 eigencoefficients
  ! yk2 - complex*16(nfreq, k) - series 2 eigencoefficients
  ! cs12 - complex*16(2*max_idx + 1, 2*max_idx + 1)
  ! is_forward - logical - if .true., calculate forward coherence
  !
    integer :: k, idx_start, idx_end, max_idx, niter, i, tmp_start, tmp_end &
      , nband, is_forward, t, neg_len
    complex*16 :: yk1(nfreq, k), yk2(nfreq, k), cs12(2*max_idx + 1, idx_end - idx_start + 1)
    complex*16, allocatable :: yk2tmp(:, :)
    ! logical :: is_forward ! change to integer for R ... ?

    niter = 2*max_idx + 1
    nband = idx_end - idx_start + 1
    allocate(yk2tmp(nband, k))

    do i = 0, niter - 1
      tmp_start = idx_start - max_idx + i
      tmp_end = tmp_start + nband - 1

!     This should take care of if f - offset occurs as a negative frequency
      if (tmp_start <= 0 .and. tmp_end <= 0) then
        yk2tmp(:, :) = conjg( yk2( (/(t, t = abs(tmp_start-2), abs(tmp_end-2), -1)/), : ) )
      else if (tmp_start <= 0 .and. tmp_end > 0) then
        neg_len = size((/(t, t = abs(tmp_start-2), 2, -1)/))
        yk2tmp(1:neg_len, :) = conjg( yk2( (/(t, t = abs(tmp_start-2), 2, -1)/), : ) )
        yk2tmp((neg_len+1):nband, :) =  yk2( 1:tmp_end, : )
      else
        yk2tmp = yk2(tmp_start:tmp_end, : )
      end if

      if (is_forward == 1) then
        cs12(i+1, :) = sum(yk1(idx_start:idx_end, : ) * conjg(yk2tmp(:, :)), 2)
      else
        cs12(i+1, :) = sum(yk1(idx_start:idx_end, : ) * yk2tmp(:, :), 2)
      end if
    end do
  end subroutine cross_spec

  subroutine block_coh(d1, d2, ndata, block_size, overlap, dt &
    , nw, k, nFFT, coh, freq, offsets, freq_range_idx &
    , max_freq_offset_idx, calc_type, is_forward)
  ! calc_type - integer - determines *what* to calculate
  !       1 - magnitude squared coherence
  !       2 - cross spectrum (complex)
  !       3 - coherency (complex valued)
  !       4 - minimum at each offset
    integer :: ndata, block_size, k, nFFT, calc_type, i &
      , freq_range_idx(2), max_freq_offset_idx, block_incr &
      , coh_nrow, coh_ncol, dstart_idx, dend_idx, is_forward &
      , ii, jj
    real*8 :: d1(ndata), d2(ndata), overlap, dt, nw &
      , df, freq(:), offsets(:)
    ! , freq_range(2), max_freq_offset not needed?
    ! freq(:) and offsets(:) were allocatable before
    real*8, allocatable :: s1(:), s2(:), s1_tot(:), s2_tot(:)
    complex*16 :: coh(:, :) ! was allocatable before
    complex*16, allocatable :: yk1(:, :), yk2(:, :), cohwrk(:, :)
    ! logical :: is_forward ! change this to an int ... for R?

    if (block_size > ndata) then
      print *, 'Block size is larger than data available.'
      return
    end if

    ! setup the fft work matrices and dpss vectors
    call fft_setup(nFFT)
    call dpss_setup(block_size, nw, k, nFFT)

    ! allocate(freq(nfreq))
    ! call freq_array(freq, nFFT, dt) ! don't need this when calling from R ?

    ! setting up for calling cross_spec() subroutine
    df = 1.0D0 / dble(dt*nFFT)
    ! max_freq_offset_idx = ceiling(max_freq_offset / df)
    ! freq_range_idx(1) = floor(freq_range(1) / df)
    ! freq_range_idx(2) = ceiling(freq_range(2) / df)

    coh_nrow = 2*max_freq_offset_idx + 1
    coh_ncol = freq_range_idx(2) - freq_range_idx(1) + 1

    ! setting up block indices
    block_incr = floor(block_size * (1.0D0 - overlap))

    ! , coh(coh_nrow, coh_ncol) was allocated
    allocate(yk1(nfreq, k), yk2(nfreq, k), cohwrk(coh_nrow, coh_ncol))
    ! , offsets(coh_nrow) was allocated
    allocate(s1(nfreq), s2(nfreq))

    if (calc_type == 2) then
      allocate(s1_tot(nfreq), s2_tot(nfreq))
      s1_tot(:) = 0.0D0
      s2_tot(:) = 0.0D0
    end if

    call offset_array(offsets, max_freq_offset_idx, df)
    i = 0
    !coh(:, :) = dcmplx(0.0D0, 0.0D0) ! was this before
    coh(:, :) = dcmplx(1.0D0, 0.0D0) ! as of Sept. 20

    do while (i*block_incr + block_size <= ndata)
      dstart_idx = i*block_incr + 1
      dend_idx = i*block_incr + block_size
      i = i + 1

      call weighted_eigencoef(d1(dstart_idx:dend_idx), block_size, dt &
       , nw, k, yk1, s1, nFFT)
      call weighted_eigencoef(d2(dstart_idx:dend_idx), block_size, dt &
       , nw, k, yk2, s2, nFFT)

      s1 = sum(abs(yk1)**2, 2)
      s2 = sum(abs(yk2)**2, 2)

      call cross_spec(yk1, yk2, cohwrk, k, freq_range_idx(1) &
       , freq_range_idx(2), max_freq_offset_idx, is_forward)

      if (calc_type == 1) then ! average MSC
        call coh_denom(cohwrk, s1, s2, coh_nrow, coh_ncol &
          , freq_range_idx(1), freq_range_idx(2), max_freq_offset_idx)
        coh = coh + dcmplx(cabs2(cohwrk, coh_nrow, coh_ncol)  , 0.0D0)
      else if (calc_type == 2) then ! average cross and auto
        coh = coh + cohwrk
        s1_tot = s1_tot + s1
        s2_tot = s2_tot + s1
      else if (calc_type == 3) then ! average coherency
        call coh_denom(cohwrk, s1, s2, coh_nrow, coh_ncol &
          , freq_range_idx(1), freq_range_idx(2), max_freq_offset_idx)
        coh = coh + cohwrk
      else if (calc_type == 4) then ! minimum MSC
        call coh_denom(cohwrk, s1, s2, coh_nrow, coh_ncol &
          , freq_range_idx(1), freq_range_idx(2), max_freq_offset_idx)
        do jj = 1, coh_ncol
          do ii = 1, coh_nrow
            !realpart(x)**2 + imagpart(x)**2
            coh(ii, jj) = dcmplx(min(realpart(coh(ii, jj))**2 + &
             imagpart(coh(ii, jj))**2 &
             , realpart(cohwrk(ii, jj))**2 + imagpart(cohwrk(ii, jj))**2) &
             , 0.0D0)
            !coh(ii, jj) = dcmplx(min(cabs2(coh(ii, jj), 1, 1), cabs2(cohwrk(ii, jj), 1, 1)), 0.0D0)
          end do
        end do
      end if
    end do

    ! don't average if we're calculating the minimum ...
    if (calc_type .ne. 4) then
      coh = coh / dble(i)
    end if

    if (calc_type == 2) then
      s1 = s1_tot / dble(i)
      s2 = s2_tot / dble(i)
      call coh_denom(coh, s1, s2, coh_nrow, coh_ncol &
        , freq_range_idx(1), freq_range_idx(2), max_freq_offset_idx)
    end if

    call fft_cleanup
    call dpss_cleanup
  end subroutine block_coh

  subroutine coh_denom(cs, s1, s2, cnrow, cncol, idx_start, idx_end, max_idx)
    integer :: idx_start, idx_end, max_idx, niter, nband &
      , tmp_start, tmp_end, cnrow, cncol, i, t, neg_len
    integer, allocatable :: idx_tmp(:)
    real*8 :: s1(nfreq), s2(nfreq)
    complex*16 :: cs(cnrow, cncol)

    niter = 2*max_idx + 1
    nband = idx_end - idx_start + 1
    allocate(idx_tmp(nband))

    do i = 0, niter - 1
      tmp_start = idx_start - max_idx + i
      tmp_end = tmp_start + nband - 1

!     This should take care of if f - offset occurs as a negative frequency
      if (tmp_start <= 0 .and. tmp_end <= 0) then
        idx_tmp = (/ (t, t = abs(tmp_start-2), abs(tmp_end-2), -1) /)
      else if (tmp_start <= 0 .and. tmp_end > 0) then
        neg_len = size((/(t, t = abs(tmp_start-2), 2, -1)/))
        idx_tmp(1:neg_len) = (/ (t, t = abs(tmp_start-2), 2, -1) /)
        idx_tmp((neg_len+1):nband) =  (/ (t, t = 1, tmp_end, 1) /)
      else
        idx_tmp = (/ (t, t = tmp_start, tmp_end, 1) /)
      end if

      cs(i+1, :) = cs(i+1, : ) / dble(sqrt((s1(idx_start:idx_end) * s2(idx_tmp))))
    end do
  end subroutine coh_denom

!! ######################################
!! ######################################
  ! subroutine tf(d1, d2, ndata, npred, block_size, overlap, dt &
  !   , nw, k, nFFT, freq_range_idx &
  !   , max_freq_offset_idx, calc_type, is_forward)
  !
  !   integer :: ndata, npred, block_size, k, nFFT, i &
  !     , freq_range_idx(2), max_freq_offset_idx(2), calc_type &
  !     , block_incr, dstart_idx, dend_idx, is_forward
  !   real*8 :: d1(ndata), d2(ndata, npred), overlap, dt, nw
  !   complex*16 :: coh(:, :) ! was allocatable before
  !   complex*16, allocatable :: yk1(:, :), yk2(:, :), cohwrk(:, :) &
  !     , H(:, :)
  !
  !   if (block_size > ndata) then
  !     print *, 'Block size is larger than data available.'
  !     return
  !   end if
  !
  !   ! setup the fft work matrices and dpss vectors
  !   call fft_setup(nFFT)
  !   call dpss_setup(block_size, nw, k, nFFT)
  !
  !   ! setting up block indices
  !   block_incr = floor(block_size * (1.0D0 - overlap))
  !
  !   allocate(yk1(nfreq, k), yk2(nfreq, k))
  !
  !
  !   call dpss_cleanup()
  !   call fft_cleanup()
  ! end subroutine tf
!! ######################################
!! ######################################

  subroutine freq_array(freq, nFFT, dt)
    integer :: i, nFFT
    real*8 :: freq(nfreq), dt

    do i = 0, nfreq - 1
      freq(i+1) = dble(i) / dble(dt*nFFT)
    end do
  end subroutine freq_array

  subroutine offset_array(offsets, max_idx, df)
    integer :: max_idx, i
    real*8 :: offsets(2*max_idx + 1), df

    do i = 1, 2*max_idx + 1
      offsets(i) = (i - max_idx - 1) * df
    end do
  end subroutine offset_array

  subroutine eigenCoef(d, ndata, nw, k, yk, nFFT, dt)
    ! Calculates the eigencoefficients and stores
    ! them in c.
    !
    ! d     - real*8    - data to be tapered and FFT'd
    ! ndata - integer*4 - number of data points
    ! nw    - real*8    - time bandwidth parameter
    ! k     - integer*4 - number of tapers
    ! c    - complex*16 - matrix to hold eigencoefficients
    ! nFFT  - integer*4 - number of frequency bins
    ! ev    - real*8    - k - eigenvalues
    !
    integer :: ndata, k, nFFT, i, col
    real*8 :: d(ndata), nw, dt
    complex*16 :: c(nFFT, k), yk(nfreq, k)
    ! logical :: isComplex

    ! isComplex = .false.

    ! calculates the slepians (double precision)
    ! call dpss(ndata, k, nw, v, ev)
    ! call dpssToEigenvalues()
    call taper_data(d, c, ndata, k, nFFT, dt)

    ! FFT each column of the tapered data matrix
    do i = 1, k
      ! print *, 'FFT loop: ', i
      call cfft1f(nFFT, inc, c(:, i), lenc, wsave, lensav, work, lenwrk, ier)
    end do

    yk = c(1:nfreq, 1:k)*nFFT
  end subroutine eigenCoef

  subroutine taper_data(d, c, n, k, nFFT, dt)
  ! Tapers the data and stores the value as a complex matrix
  ! Ready to be FFT'd
  !
  ! d   - real*8(n)  - data to be tapered
  ! c   - complex*16(nFFT, k) - matrix to hold tapered data (real part)
  ! n   - integer   - length of data vector
  ! k   - integer   - number of data tapers to use
  ! nFFT  - integer - number of frequency bins to use (zeropad)
  !
    integer :: n, k, nFFT, i, j
    real*8 :: d(n), dt
    complex*16 :: c(nFFT, k)

    ! could change this to use:
    ! c(:, :) = dcmplx(0.0D0, 0.0D0)
    do j = 1, k
      c(1:n, j) = d * (sqrt(dt) * v(:, j))
      do i = n+1, nFFT
        c(i, j) = dcmplx(0.0D0, 0.0D0)
      end do
    end do
  end subroutine taper_data

  subroutine dpssToEigenvalues(ndata, k, nw, nFFT)
  !+++++++++++++++++++++++++++++++++
  ! Calculates eigenvalues that will be between 0 and 1
  ! Modified from multitaper function dpssToEigenvalues()
  ! Percival and Walden (1993) exercise 8.5
  !
  ! ndata - integer - length of the tapers
  ! k     - integer - number of tapers
  ! nw    - real*8  - time bandwidth parameter
  ! nFFT  - integer - number of frequency bins to use (zeropadding)
  !
    integer :: ndata, k, nFFT, i, j, l !, npot
    real*8 :: nw, w, ratios(ndata), evtmp
    complex*16 :: cx(nFFT, k)

    w = nw / dble(ndata)
    ! use nFFT ?  npot = 2**ceiling(log(dble(2*ndata))/log(2.0D0))
    cx(:, :) = dcmplx(0.0D0, 0.0D0)
    cx(1:ndata, :) = dcmplx(v(:, :))

    ratios(1) = 2.0D0*w
    do j = 1, ndata-1
      ratios(j+1) = sin(2*pi*w*j)/(pi*j)
    end do

    do j = 1, k
      call cfft1f(nFFT, inc, cx(:, j), lenc, wsave, lensav, work, lenwrk, ier)

      cx(:, j) = dcmplx(abs(cx(:, j)*nFFT)**2, 0.0D0)

      call cfft1f(nFFT, inc, cx(:, j), lenc, wsave, lensav, work, lenwrk, ier)

      ev(j) = sum(realpart(cx(2:ndata, j)) * (ratios(2:ndata)))
      ! evtmp = 0.0D0
      ! do i = ndata, 2, -1
      !   evtmp = evtmp + realpart(cx(i, j)) * ratios(i)
      ! end do
      !
      ! ev(j) = 2*evtmp + ratios(1)*realpart(cx(1, j))
      ev(j) = 2*ev(j) + ratios(1)*realpart(cx(1, j))
    end do
  end subroutine dpssToEigenvalues

  subroutine adaptive_weights(yk, dk, k, spec, dofs, dofav, var, dt)
  ! implement a wrapper for mw2wta to get the weights only
  !
  ! yk  - complex*16(nfreq, k)  - eigencoefficients
  ! dk  - real*8(nfreq, k)      - matrix to hold weights
  ! k   - integer               - number of tapers
  ! spec  - real*8(nfreq)       - multitaper spectral estimate
  ! dofs  - real*8(nfreq)       - degrees of freedom at each frequency
  ! dofav - real*8              - average dofs across all freqs?
  ! var   - real*8              - variance of the time series
  ! dt    - real*8              - sampling rate
  !
    integer :: k, maxadit, mxiter
    real*8 :: sa(nfreq, k), dofs(nfreq)
    real*8 :: dofav, dk(nfreq, k), spec(nfreq), evp(k), tol, aviter, var, dt
    complex*16 :: yk(nfreq, k)

    sa = dble(abs(yk)**2)
    evp = 1.0D0 - ev
    tol = 0.03D0
    maxadit = 100

    call mw2wta(sa, dk, nfreq, k, spec, ev, evp, dofs, dofav, var, dt, tol, maxadit, mxiter, aviter)

    dk = dsqrt(dk)
  end subroutine adaptive_weights

  subroutine variance(d, ndata, var)
  ! calculates the standard unbiased variance
  !
  ! d - real*8(ndata) - time series
  ! ndata - integer   - length of the time series
  ! var   - real*8    - variance of the time series
  !
    integer :: ndata
    real*8 :: d(ndata), var, mean

    mean = sum(d) / dble(ndata)
    var = 0.0D0
    ! biased sample variance!! (copying pgk:multitaper)
    var = sum((d - mean)**2) / (ndata)
  end subroutine variance

  subroutine weighted_eigencoef(d, ndata, dt, nw, k, yk, spec, nFFT)
  ! calculates and returns the weighted eigencoefficients
  !
  ! d   - real*8(ndata) - time series (data)
  ! ndata - integer     - length of the time series data
  ! dt    - real*8      - sampling rate
  ! nw    - real*8      - time bandwidth parameter
  ! k     - integer     - number of tapers to use
  ! yk - complex*16(nfreq, k) - matrix to hold the weighted eigencoefficients
  ! nFFT  - integer     - number of frequency bins to use (zeropadding)
  !
    integer :: ndata, nFFT, k
    real*8 :: d(ndata), nw, dt, var, dk(nfreq, k), spec(nfreq), dofs(nfreq), dofav
    complex*16 :: yk(nfreq, k)

    call variance(d, ndata, var)
    call eigencoef(d, ndata, nw, k, yk, nFFT, dt)
    call adaptive_weights(yk, dk, k, spec, dofs, dofav, var, dt)

    yk = yk * dk
  end subroutine weighted_eigencoef

  subroutine zSvdRegression(Y, X, m, n, beta, stdErr, svd_ev)
  ! uses the SVD of the matrix X to calculate the regression coefficients for
  ! the model: Y = XB (where B is a vector of beta coefficients)
  !
  ! as a side note: this mxn notation is bad... as a result of lapack docs
  !
  ! [in] Y      - complex*16(m)  - response vector
  ! [in] X      - complex*16(m,n)  - design matrix
  ! [in] m      - integer - number of rows of both Y and X
  ! [in] n      - integer - number of columns of X
  ! [out] beta  - complex*16(n) - the regression coefficients
  ! [out] stdErr  - real*8(n) - the standard error on the beta estimates
  ! [out] ev    - real*8(n) - the square of the singular values

    integer :: i, m, n, lda, ldu, ldvt, lwork, info, lrwork
    real*8 :: svd_ev(n), stdErr(n), lworkopt
    real*8, allocatable :: rwork(:), s(:)
    complex*16 :: Y(m), X(m, n), beta(n)
    complex*16, allocatable :: svd_work(:), u(:, :), vt(:, :)
    character(1) :: jobu, jobvt

    ! 'S' says that we want the left and right singular vectors
    jobu = 'S'
    jobvt = 'S'
    lda = m
    ldu = m
    ldvt = n

    ! set values according to:
    ! http://www.netlib.org/lapack/explore-html/index.html
    ! search for zgesvd
    lrwork = 5*min(m, n)
    allocate(rwork(lrwork))
    allocate(s(min(m,n)))
    allocate(u(ldu, min(m,n)))
    allocate(vt(ldvt, n))

    ! obtain optimal size for lwork
    lwork = -1
    call zgesvd(jobu, jobvt, m, n, X, lda, s, u, ldu, vt, ldvt &
      , lworkopt, lwork, rwork, info)

    ! allocate the work array
    lwork = nint(lworkopt)
    allocate(svd_work(lwork))

    ! perform the svd
    call zgesvd(jobu, jobvt, m, n, X, lda, s, u, ldu, vt, ldvt &
      , svd_work, lwork, rwork, info)

    ! calculate them betas '*' is matrix mult here
    ! beta = V * s^(-1) * [ t(u) Y ] (mandel - eqn's (17) and (12))
    beta = matmul( transpose(conjg(vt)), matmul( transpose(conjg(u)), Y ) / s )

    ! calculate the standard error estiamtes on the betas
    ! mandel eqn (20) NOT YET IMPLEMENTED
    do i = 1, n
      stdErr(i) = -1 !sum(vt(:, i)x)
    end do

    do i = 1, n
      svd_ev(i) =  -1 ! NOT YET IMPLEMENTED!
    end do
  end subroutine zSvdRegression

  function cabs2(x, nrow, ncol)
    integer :: nrow, ncol
    real*8 :: cabs2(nrow, ncol)
    complex*16 :: x(nrow, ncol)

    cabs2 = realpart(x)**2 + imagpart(x)**2
  end function cabs2

  subroutine testem(x)
    integer :: x

    x = 5
  end subroutine testem
end module mtm_mod
