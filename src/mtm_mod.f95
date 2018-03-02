module mtm_mod
  implicit none
  real*8, parameter :: pi = dble(3.1415926535897932)
  real*8, allocatable :: wsave(:), work(:), v(:, :), ev(:), &
    wsave2(:), work2(:), v2(:, :), ev2(:)
  integer :: lenc, lensav, lenwrk, inc, ier, nfreq, &
    lenc2, lensav2, lenwrk2, inc2, ier2, nfreq2
  logical :: is_dpss_setup = .false., is_fft_setup = .false.

contains
  subroutine fft_setup(nFFT, id)
  ! calls the setup routine for fftpack5.1d routines
  ! cfft1f (or cfft1b)
  ! id - integer - indicates whether data set 1 or 2 due to
  ! different sampling rates!
    integer :: nFFT, id

    if (id .eq. 1) then
      lenc = nFFT
      nfreq = nFFT/2 + 1
      inc = 1
      lensav = 2*nFFT + int(log( dble(nFFT) ) / log( dble(2) )) + 4
      lenwrk = 2*nFFT

      allocate(wsave(lensav))
      allocate(work(lenwrk))

      call cfft1i(nFFT, wsave, lensav, ier)
    else if (id .eq. 2) then
      lenc2 = nFFT
      nfreq2 = nFFT/2 + 1
      inc2 = 1
      lensav2 = 2*nFFT + int(log( dble(nFFT) ) / log( dble(2) )) + 4
      lenwrk2 = 2*nFFT

      allocate(wsave2(lensav2))
      allocate(work2(lenwrk2))

      call cfft1i(nFFT, wsave2, lensav2, ier2)
    end if
  end subroutine fft_setup

  subroutine fft_cleanup(id)
  ! Frees up memory again so that this code can be re-called.
    integer :: id

    if (id .eq. 1) then
      deallocate(wsave, work)
    else if (id .eq. 2) then
      deallocate(wsave2, work2)
    end if
  end subroutine

  subroutine dpss_setup(ndata, nw, k, nFFT, id)
  ! calculates the dpss tapers and eigenvalues
  !
  ! ndata - integer - length of data to be used
  ! nw    - real*8  - time bandwidth parameter
  ! k     - integer - number of tapers to be used
  ! nFFT  - integer - length to zeropad out to
  ! id    - integer - which dataset to use (currently 1 or 2)
  !
    integer :: ndata, k, nFFT, id
    real*8 :: nw

    if (id .eq. 1) then
      allocate(v(ndata, k), ev(k))
      call dpss(ndata, k, nw, v, ev)
    else if (id .eq. 2) then
      allocate(v2(ndata, k), ev2(k))
      call dpss(ndata, k, nw, v2, ev2)
    end if

    call dpssToEigenvalues(ndata, k, nw, nFFT, id)
  end subroutine dpss_setup

  subroutine dpss_cleanup(id)
  ! Frees up memory again so that this code can be re-called.
    integer :: id

    if (id .eq. 1) then
      deallocate(v, ev)
    else if (id .eq. 2) then
      deallocate(v2, ev2)
    end if
  end subroutine dpss_cleanup

  subroutine cross_spec(yk1, yk2, cs12, k, idx_start, idx_end, max_idx &
    , fRatio, is_forward)
  ! Calculates the cross spectrum given two eigencoefficient matrices
  !
  ! yk1 - complex*16(nfreq, k) - series 1 eigencoefficients
  ! yk2 - complex*16(nfreq2, k) - series 2 eigencoefficients
  ! cs12 - complex*16(2*max_idx + 1, 2*max_idx + 1)
  ! is_forward - logical - if .true., calculate forward coherence
  !
    integer :: k, idx_start, idx_end, max_idx, niter, i, tmp_start, tmp_end &
      , nband, is_forward, t, neg_len, fRatio
    complex*16 :: yk1(nfreq, k), yk2(nfreq2, k) &
      , cs12(2*max_idx + 1, idx_end - idx_start + 1)
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
        yk2tmp(:, :) = conjg( yk2( (/(t, t = fRatio * abs(tmp_start-2) &
          - (fRatio - 1) &
          , fRatio * abs(tmp_end-2) - (fRatio - 1), -1*fRatio)/), : ) )
      else if (tmp_start <= 0 .and. tmp_end > 0) then
        neg_len = size((/(t, t = abs(tmp_start-2), 2, -1)/))
        yk2tmp(1:neg_len, :) = conjg( yk2( (/(t, t = fRatio * abs(tmp_start-2) &
          - (fRatio-1), fRatio + 1, -1*fRatio)/), : ) )
        yk2tmp((neg_len+1):nband, :) =  yk2( (/(t, t = 1 &
          , fRatio*tmp_end - (fRatio - 1))/), : )
      else
        yk2tmp = yk2( (/ (t, t = fRatio*tmp_start - (fRatio - 1) &
          , fRatio*tmp_end - (fRatio - 1), fRatio) /), : )
      end if

      if (is_forward == 1) then
        cs12(i+1, :) = sum(yk1(idx_start:idx_end, : ) * conjg(yk2tmp(:, :)), 2)
      else
        cs12(i+1, :) = sum(yk1(idx_start:idx_end, : ) * yk2tmp(:, :), 2)
      end if
    end do
  end subroutine cross_spec

  subroutine block_coh(d1, d2, ndata, ndata2, block_size, block_size2 &
    , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio, coh, freq &
    , offsets, freq_range_idx, max_freq_offset_idx &
    , calc_type, conv_msc2norm, is_forward)
  ! calc_type - integer - determines *what* to calculate
  !       1 - magnitude squared coherence
  !       2 - cross spectrum (complex)
  !       3 - coherency (complex valued)
  !       4 - minimum at each offset
  !
  ! conv_msc2norm - integer - whether to convert the MSC to Standard normal
  !     before averaging or what-have you - ONLY if calc_type == 1
    integer :: ndata, ndata2, block_size, block_size2 &
      , k, nFFT, nFFT2, calc_type, i &
      , freq_range_idx(2), max_freq_offset_idx &
      , block_incr, block_incr2 &
      , coh_nrow, coh_ncol, dstart_idx, dstart_idx2 &
      , dend_idx, dend_idx2, conv_msc2norm, is_forward &
      , ii, jj, fRatio, dtRatio
    real*8 :: d1(ndata), d2(ndata2), overlap, dt, dt2, nw, nw2 &
      , df, df2, freq(:), offsets(:)
    ! , freq_range(2), max_freq_offset not needed?
    ! freq(:) and offsets(:) were allocatable before
    real*8, allocatable :: s1(:), s2(:), s1_tot(:), s2_tot(:), mscwrk(:, :)
    complex*16 :: coh(:, :) ! was allocatable before
    complex*16, allocatable :: yk1(:, :), yk2(:, :), cohwrk(:, :)
    ! logical :: is_forward ! change this to an int ... for R?

    if (block_size > ndata) then
      print *, 'Block size is larger than data available.'
      return
    end if

    ! setup the fft work matrices and dpss vectors
    call fft_setup(nFFT, 1)
    call fft_setup(nFFT2, 2)
    call dpss_setup(block_size, nw, k, nFFT, 1)
    call dpss_setup(block_size2, nw2, k, nFFT2, 2)

    ! fRatio = (nFFT2 * dt2) / (nFFT * dt)
    dtRatio = dt / dt2

    ! allocate(freq(nfreq))
    ! call freq_array(freq, nFFT, dt) ! don't need this when calling from R ?

    ! setting up for calling cross_spec() subroutine
    df = 1.0D0 / dble(dt*nFFT)
    df2 = 1.0D0 / dble(dt2*nFFT2)
    ! max_freq_offset_idx = ceiling(max_freq_offset / df)
    ! freq_range_idx(1) = floor(freq_range(1) / df)
    ! freq_range_idx(2) = ceiling(freq_range(2) / df)

    coh_nrow = 2*max_freq_offset_idx + 1
    coh_ncol = freq_range_idx(2) - freq_range_idx(1) + 1

    ! setting up block indices
    block_incr = floor(block_size * (1.0D0 - overlap))
    block_incr2 = block_incr * dtRatio

    ! , coh(coh_nrow, coh_ncol) was allocated
    allocate(yk1(nfreq, k), yk2(nfreq2, k), cohwrk(coh_nrow, coh_ncol))
    ! , offsets(coh_nrow) was allocated
    allocate(s1(nfreq), s2(nfreq2))

    if (calc_type .eq. 1) then
      allocate(mscwrk(coh_nrow, coh_ncol))
    end if

    if (calc_type == 2) then
      allocate(s1_tot(nfreq), s2_tot(nfreq2))
      s1_tot(:) = 0.0D0
      s2_tot(:) = 0.0D0
    end if

    call offset_array(offsets, max_freq_offset_idx, df)
    i = 0
    if (calc_type == 4) then
      coh(:, :) = dcmplx(1.0D0, 0.0D0) ! as of Sept. 20
    else
      coh(:, :) = dcmplx(0.0D0, 0.0D0) ! was this before
    end if

    do while (i*block_incr + block_size <= ndata)
      dstart_idx = i*block_incr + 1
      dend_idx = i*block_incr + block_size - 1
      dstart_idx2 = i*block_incr2 + 1
      dend_idx2 = i*block_incr2 + block_size2 - 1
      i = i + 1

      call weighted_eigencoef(d1(dstart_idx:dend_idx), block_size, dt &
       , nw, k, yk1, s1, nFFT, 1)
      call weighted_eigencoef(d2(dstart_idx2:dend_idx2), block_size2, dt2 &
       , nw2, k, yk2, s2, nFFT2, 2)

      s1 = sum(abs(yk1)**2, 2)
      s2 = sum(abs(yk2)**2, 2) ! multiply by fRatio here - maybe?

      call cross_spec(yk1, yk2, cohwrk, k, freq_range_idx(1) &
       , freq_range_idx(2), max_freq_offset_idx, fRatio, is_forward)

      if (calc_type == 1) then ! average MSC
        call coh_denom(cohwrk, s1, s2, coh_nrow, coh_ncol &
          , freq_range_idx(1), freq_range_idx(2), max_freq_offset_idx, fRatio)
        if (conv_msc2norm .eq. 1) then
          mscwrk = cabs2(cohwrk, coh_nrow, coh_ncol)
          call msc2norm(mscwrk, coh_nrow, coh_ncol, k)
          coh = coh + dcmplx( mscwrk , 0.0D0 )
        else
          coh = coh + dcmplx(cabs2(cohwrk, coh_nrow, coh_ncol)  , 0.0D0)
        end if
      else if (calc_type == 2) then ! average cross and auto
        coh = coh + cohwrk
        s1_tot = s1_tot + s1
        s2_tot = s2_tot + s1
      else if (calc_type == 3) then ! average coherency
        call coh_denom(cohwrk, s1, s2, coh_nrow, coh_ncol &
          , freq_range_idx(1), freq_range_idx(2), max_freq_offset_idx, fRatio)
        coh = coh + cohwrk
      else if (calc_type == 4) then ! minimum MSC
        call coh_denom(cohwrk, s1, s2, coh_nrow, coh_ncol &
          , freq_range_idx(1), freq_range_idx(2), max_freq_offset_idx, fRatio)
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
      else
        coh = dcmplx(999.0D0, 999.0D0)
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
        , freq_range_idx(1), freq_range_idx(2), max_freq_offset_idx, fRatio)
    end if

    call fft_cleanup(1)
    call fft_cleanup(2)
    call dpss_cleanup(1)
    call dpss_cleanup(2)
  end subroutine block_coh

  subroutine coh_denom(cs, s1, s2, cnrow, cncol, idx_start, idx_end &
    , max_idx, fRatio)
    integer :: idx_start, idx_end, max_idx, niter, nband &
      , tmp_start, tmp_end, cnrow, cncol, i, t, neg_len &
      , fRatio
    integer, allocatable :: idx_tmp(:)
    real*8 :: s1(nfreq), s2(nfreq2)
    complex*16 :: cs(cnrow, cncol)

    niter = 2*max_idx + 1
    nband = idx_end - idx_start + 1
    allocate(idx_tmp(nband))

    do i = 0, niter - 1
      tmp_start = idx_start - max_idx + i
      tmp_end = tmp_start + nband - 1

!     This should take care of if f - offset occurs as a negative frequency
      if (tmp_start <= 0 .and. tmp_end <= 0) then
        idx_tmp = (/ (t, t = fRatio * abs(tmp_start-2) - (fRatio - 1) &
          , fRatio * abs(tmp_end-2) - (fRatio - 1), -1*fRatio) /)
      else if (tmp_start <= 0 .and. tmp_end > 0) then
        neg_len = size((/(t, t = abs(tmp_start-2), 2, -1)/))
        idx_tmp(1:neg_len) = (/ (t &
          , t = fRatio * abs(tmp_start-2) - (fRatio - 1), fRatio + 1 &
          , -1*fRatio) /)
        idx_tmp((neg_len+1):nband) =  (/ (t, t = 1 &
          , fRatio * tmp_end - (fRatio - 1), 1*fRatio) /)
      else
        idx_tmp = (/ (t, t = fRatio * tmp_start - (fRatio - 1) &
          , fRatio * tmp_end - (fRatio - 1), 1*fRatio) /)
      end if

      cs(i+1, :) = cs(i+1, : ) / dble(sqrt((s1(idx_start:idx_end) * s2(idx_tmp))))
    end do
  end subroutine coh_denom

  subroutine tf_zero(d1, d2, ndata, ndata2, npred, block_size, block_size2 &
    , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio &
    , freq_range_idx, H, n_row_H)
  ! [in] d1(ndata) - real*8 - a vector of data (response)
  ! [in] d2(ndata2, npred) - real*8 - predictors - each column is one series
  ! [in] ndata - integer - number of total data points in response
  ! [in] ndata2 - integer - number of total data points in each predictor
  ! [in] npred - integer - number of predictors
  ! [in] block_size - integer - block size for series1
  ! [in] block_size2 - integer - block size for series2
  ! [in] overlap - real*8 - proportion that each block overlaps
  ! [in] dt - real*8 - sampling rate for series1
  ! [in] dt2 - real*8 - sampling rate for series2
  ! [in] nw - real*8 - time bandwidth parameter for multitaper, series1
  ! [in] nw2 - real*8 - time bandwidth parameter for multitaper, series2
  ! [in] k - integer - number of tapers (2*NW-1 usually)
  ! [in] nFFT - integer - number of frequency bins (zero-padding), series1
  ! [in] nFFT2 - integer - number of frequency bins (zero-padding), series2
  ! [in] fRatio - integer - spacing between series 2 frequency bins such that
  ! subsequent frequencies for series1 are the same as frequencies separated by
  ! fRatio in series2
  ! [in] freq_range_idx - integer(2) - starting and ending indices over which
  ! to estimate the transfer function
  ! [out] H - complex*16(n_row_H, npred) - transfer function for each predictor
  ! in each column.

    integer :: ndata, ndata2, npred, block_size, block_size2, k &
      , nFFT, nFFT2, fRatio, block_incr, block_incr2, nblocks, nblocks2 &
      , dtRatio, freq_range_idx(2), n_row_H &
      , i, j, p
    real*8 :: d1(ndata), d2(ndata2, npred), overlap, dt, dt2, nw, nw2 &
      , stdErr(npred), eigenval(npred)
    complex*16 :: H(n_row_H, npred)
    complex*16, allocatable :: yk1(:, :, :), yk2(:, :, :, :) &
      , Y(:), design(:, :)

    dtRatio = dt / dt2 ! doesn't check that this should *actually* be an integer

    block_incr = block_increment(block_size, overlap)
    block_incr2 = block_incr * dtRatio

    nblocks = calculate_nblocks(ndata, block_size, block_incr)

    call fft_setup(nFFT, 1)
    call fft_setup(nFFT2, 2)
    call dpss_setup(block_size, nw, k, nFFT, 1)
    call dpss_setup(block_size2, nw2, k, nFFT2, 2)

    allocate(design(nblocks*k, npred), Y(nblocks*k))
    allocate(yk1(nfreq, k, nblocks), yk2(nfreq2, k, npred, nblocks))

    call calc_tf_wt_eigen(block_incr, block_incr2 &
      , block_size, block_size2, nblocks, d1, d2, npred &
      , dt, dt2, nw, nw2, k, nFFT, nFFT2, yk1, yk2)

    do j = freq_range_idx(1), freq_range_idx(2)
      do i = 0, nblocks-1
        Y((i*k+1):(i+1)*k) = yk1(j, :, i+1)
        do p = 1, npred
          design((i*k+1):(i+1)*k, p) = yk2(j, :, p, i+1)
        end do
      end do
      call zSvdRegression(Y, design, nblocks*k, npred &
        ,H(j, :), stdErr, eigenval)
    end do

    call dpss_cleanup(1)
    call dpss_cleanup(2)
    call fft_cleanup(1)
    call fft_cleanup(2)
  end subroutine tf_zero

!! ######################################
!! ######################################
  subroutine tf(d1, d2, ndata, ndata2, npred, block_size, block_size2 &
    , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio &
    , freq_range_idx, max_freq_offset_idx, H &
    , coh_nrow, coh_ncol, totFreqByCol, totFreqByColByPred, total_offsets &
    , col_order, hPredBreak, hIdx, nhIdx)

    integer :: ndata, ndata2, npred, nblocks, nblocks2 &
      , block_size, block_size2, k, nFFT, nFFT2 &
      , freq_range_idx(2), max_freq_offset_idx &
      , block_incr, block_incr2, dstart_idx, dend_idx &
      , dstart_idx2, dend_idx2  &
      , fRatio, dtRatio, use_off, nOff &
      , total_offsets &
      , coh_nrow, coh_ncol &
      , col_order(coh_nrow, coh_ncol, npred) &
      , totFreqByCol(coh_ncol) &
      , offIdxAll(coh_nrow) &
      , idxAll(coh_nrow), negIdxAll(3), nNegOff, nPosOff &
      , i, j, p, ii, desCol(npred), t, wrk2(coh_nrow) &
      , hIdx(nhIdx), nhIdx, nhSubIdx, hPredBreak(npred+1) &
      , totFreqByColByPred(npred+1, coh_ncol)
    integer, allocatable :: freqIdx2(:), negOffIdx(:), posOffIdx(:) &
      , hSubIdx(:)
    real*8 :: d1(ndata), d2(ndata2, npred), overlap &
      , dt, dt2, nw, nw2, freq(nFFT/2+1)
    real*8, allocatable :: tmpErr(:), tmpSvd_Ev(:)
    complex*16 :: H(coh_ncol, total_offsets)
    complex*16, allocatable :: yk1(:, :, :), yk2(:, :, :, :) &
      , design(:, :), Y(:), Htmp(:)

    ! fRatio = (nFFT2 *dt2) / (nFFT * dt) <- gets passed
    dtRatio = dt / dt2 ! I dont think this gets used...

    ! setting up block indices for eigencoefficient calcs
    block_incr = block_increment(block_size, overlap) !floor(block_size * (1.0D0 - overlap))
    block_incr2 = block_incr * dtRatio
  ! these two sizes should be the same...
    ! nblocks = min(size( (/ (j, j = 1, ndata - block_size + 1, block_incr) /) ) &
    !   , size( (/ (j, j = 1, ndata2 - block_size2 + 1, block_incr2) /) ))
    nblocks = calculate_nblocks(ndata, block_size, block_incr)
    nblocks2 = calculate_nblocks(ndata2, block_size2, block_incr2)

!!! make this work with R!
    if (nblocks .ne. nblocks2) then
      print *, "I don't know what to do here... error number?"
    end if

  ! setup the fft work matrices and dpss vectors
    call fft_setup(nFFT, 1)
    call fft_setup(nFFT2, 2)
    call dpss_setup(block_size, nw, k, nFFT, 1)
    call dpss_setup(block_size2, nw2, k, nFFT2, 2)

    allocate(yk1(nfreq, k, nblocks), yk2(nfreq2, k, npred, nblocks))
    allocate(freqIdx2(nfreq2))

    call calc_tf_wt_eigen(block_incr, block_incr2 &
      , block_size, block_size2, nblocks, d1, d2, npred &
      , dt, dt2, nw, nw2, k, nFFT, nFFT2, yk1, yk2)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! setup the design matrices and allocate H(:, :)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! idxAll = (/ (i, i = 1, coh_nrow) /) ! <- delete this
    idxAll = (/ (i, i = freq_range_idx(1), freq_range_idx(2)) /)

    !allocate(H(size(idxAll), total_offsets)) ! <- gets passed in
    H = dcmplx(0.0D0, 0.0D0) ! initialize this
    allocate(Y(nblocks*k))

    do j = 1, coh_ncol
      ! no significant frequencies to use at this central frequency
      if (totFreqByCol(j) == 0) then
        cycle
      end if

      ! current offset indices in yk2 to be used
      ! central frequency - offset up to central frequency + offset
      offIdxAll = (/ (ii, ii = idxAll(j) - max_freq_offset_idx &
        , idxAll(j) + max_freq_offset_idx) /)

      ! holds the indices of columns of H for this central frequency that
      ! will contain values
      allocate(hSubIdx(totFreqByCol(j)), Htmp(totFreqByCol(j)))

      ! design matrix changes size based on number of significant offsets
      allocate( design( nblocks*k, totFreqByCol(j) ) )

      ! vector containing start and end indices for putting
      ! covariates from each predictor
      desCol = (/ 1, (0, t = 1, npred) /)
      do p = 1, npred

        ! don't waste time if there are no frequencies to add from this
        ! predictor
        if (totFreqByColByPred(p+1, j) == 0) then
          cycle
        end if
        ! which offsets for the j-th column and p-th predictor should be used
        wrk2 = col_order(:, j, p)
        ! the number of predictors to assign indices to the design matrix
        ! desCol(p+1) = desCol(p) + size(pack(offIdxAll, wrk2 > 0)) <- wrong?
        desCol(p+1) = desCol(p) + size(pack(wrk2, wrk2 > 0)) ! new

      ! determine number of offsets that are at negative frequencies
      ! and therefore require conjugation :: X(-f) = X*(f)
        nNegOff = size(pack(offIdxAll, (wrk2 > 0) .and. (offIdxAll .le. 0) ))
        nPosOff = size(pack(offIdxAll, (wrk2 > 0) .and. (offIdxAll > 0) ))
        allocate(negOffIdx(nNegOff), posOffIdx(nPosOff))
        ! takes into account fRatio :)
        negOffIdx = ( fRatio * &
          (pack(offIdxAll, (wrk2 > 0) .and. (offIdxAll .le. 0) ) - 1) ) + 1
        posOffIdx = (fRatio * &
          (pack(offIdxAll, (wrk2 > 0) .and. (offIdxAll > 0) ) - 1) ) + 1

        ! these are columns of the design matrix which will have the -ve offs
        negIdxAll = (/ 1, nNegOff, desCol(p+1) - desCol(p) /)

        do i = 0, nblocks-1
          ! if there *are* offsets at negative frequencies:
          if (negIdxAll(2) .ne. 0) then
            design((i*k+1):((i+1)*k), desCol(p):(desCol(p) + nNegOff - 1)) = &
              conjg(transpose(yk2(negOffIdx, :, p, i+1)))
          end if
          ! if there are positive offset frequencies
          if (negIdxAll(3) .ne. negIdxAll(2)) then
            design((i*k+1):((i+1)*k), (desCol(p)+nNegOff):(desCol(p+1)-1)) = &
              transpose(yk2(posOffIdx, :, p, i+1))
          end if
          if (i .eq. 0) then
            ! response eigencoefficients
            Y((i*k+1):((i+1)*k)) = yk1(freq_range_idx(1) + j - 1, :, i+1)
          end if

          deallocate(negOffIdx, posOffIdx)
        end do

        ! determind which columns of H should be used for this central freq
        call findHidx(col_order(:, j, p) &
          , hIdx(hPredBreak(p):(hPredBreak(p+1)-1)) &
          , coh_nrow & !, size(hPredBreak(p):(hPredBreak(p+1)-1)) &
          , hSubIdx((sum(totFreqByColByPred(1:p &
            , j))+1):(sum(totFreqByColByPred(1:(p+1), j)))) &
          , hPredBreak(p)-1) !, totFreqByColByPred(p, j))
      end do

      call zSvdRegression(Y, design, nblocks*k, totFreqByCol(j), Htmp &
        , tmpErr, tmpSvd_Ev)

      H(j, hSubIdx) = Htmp

      deallocate(design, hSubIdx, Htmp)
    end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cleanup
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call dpss_cleanup(1)
    call dpss_cleanup(2)
    call fft_cleanup(1)
    call fft_cleanup(2)
  end subroutine tf
!! ######################################
!! ######################################

!!! This doesn't ever get called - did not change to work with differing
! sample rates between series.
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

  subroutine eigenCoef(d, ndata, nw, k, yk, nFFT, dt, id)
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
    integer :: ndata, k, nFFT, i, col, id
    real*8 :: d(ndata), nw, dt
    complex*16 :: c(nFFT, k)
    complex*16 :: yk(:, :)

    if (id .eq. 1) then
      call taper_data(d, c, ndata, k, nFFT, dt, id)

      ! FFT each column of the tapered data matrix
      do i = 1, k
        ! print *, 'FFT loop: ', i
        call cfft1f(nFFT, inc, c(:, i), lenc, wsave, lensav, work, lenwrk, ier)
      end do

      yk = c(1:nfreq, 1:k)*nFFT
  ! second series
    else if (id .eq. 2) then
      call taper_data(d, c, ndata, k, nFFT, dt, id)

      ! FFT each column of the tapered data matrix
      do i = 1, k
        ! print *, 'FFT loop: ', i
        call cfft1f(nFFT, inc2, c(:, i), lenc2, wsave2, lensav2 &
          , work2, lenwrk2, ier2)
      end do

      yk = c(1:nfreq2, 1:k)*nFFT
    end if
  end subroutine eigenCoef

  subroutine taper_data(d, c, n, k, nFFT, dt, id)
  ! Tapers the data and stores the value as a complex matrix
  ! Ready to be FFT'd
  !
  ! d   - real*8(n)  - data to be tapered
  ! c   - complex*16(nFFT, k) - matrix to hold tapered data (real part)
  ! n   - integer   - length of data vector
  ! k   - integer   - number of data tapers to use
  ! nFFT  - integer - number of frequency bins to use (zeropad)
  !
    integer :: n, k, nFFT, i, j, id
    real*8 :: d(n), dt
    complex*16 :: c(nFFT, k)

    ! could change this to use:
    ! c(:, :) = dcmplx(0.0D0, 0.0D0)
  ! first series
    if (id .eq. 1) then
      do j = 1, k
        c(1:n, j) = d * (sqrt(dt) * v(:, j))
        do i = n+1, nFFT
          c(i, j) = dcmplx(0.0D0, 0.0D0)
        end do
      end do
  ! second series
    else if (id .eq. 2) then
      do j = 1, k
        c(1:n, j) = d * (sqrt(dt) * v2(:, j))
        do i = n+1, nFFT
          c(i, j) = dcmplx(0.0D0, 0.0D0)
        end do
      end do
    end if
  end subroutine taper_data

  subroutine dpssToEigenvalues(ndata, k, nw, nFFT, id)
  !+++++++++++++++++++++++++++++++++
  ! Calculates eigenvalues that will be between 0 and 1
  ! Modified from multitaper function dpssToEigenvalues()
  ! Percival and Walden (1993) exercise 8.5
  !
  ! ndata - integer - length of the tapers
  ! k     - integer - number of tapers
  ! nw    - real*8  - time bandwidth parameter
  ! nFFT  - integer - number of frequency bins to use (zeropadding)
  ! id    - integer - which dataset to use (currently 1 or 2)
  !
    integer :: ndata, k, nFFT, i, j, l, id !, npot
    real*8 :: nw, w, ratios(ndata), evtmp
    complex*16 :: cx(nFFT, k)

    w = nw / dble(ndata)
    ! use nFFT ?  npot = 2**ceiling(log(dble(2*ndata))/log(2.0D0))
    cx(:, :) = dcmplx(0.0D0, 0.0D0)

    ratios(1) = 2.0D0*w
    do j = 1, ndata-1
      ratios(j+1) = sin(2*pi*w*j)/(pi*j)
    end do

    if (id .eq. 1) then
      cx(1:ndata, :) = dcmplx(v(:, :))

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
    else if (id .eq. 2) then
      cx(1:ndata, :) = dcmplx(v2(:, :))

      do j = 1, k
        call cfft1f(nFFT, inc2, cx(:, j), lenc2, wsave2, lensav2, work2, lenwrk2, ier2)

        cx(:, j) = dcmplx(abs(cx(:, j)*nFFT)**2, 0.0D0)

        call cfft1f(nFFT, inc2, cx(:, j), lenc2, wsave2, lensav2, work2, lenwrk2, ier2)

        ev2(j) = sum(realpart(cx(2:ndata, j)) * (ratios(2:ndata)))
        ev2(j) = 2*ev2(j) + ratios(1)*realpart(cx(1, j))
      end do
    end if
  end subroutine dpssToEigenvalues

  subroutine adaptive_weights(yk, dk, k, spec, dofs, dofav, var, dt, id)
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
    integer :: k, maxadit, mxiter, id
    real*8, allocatable :: sa(:, :)
    real*8 :: dofav, evp(k), tol, aviter, var, dt, spec(:), dk(:, :), dofs(:)
    complex*16 :: yk(:, :)

    if (id .eq. 1) then
      allocate(sa(nfreq, k))
      sa = dble(abs(yk)**2)
      evp = 1.0D0 - ev
      tol = 0.03D0
      maxadit = 100
      call mw2wta(sa, dk, nfreq, k, spec, ev, evp, dofs, dofav, var &
        , dt, tol, maxadit, mxiter, aviter)
    else if (id .eq. 2) then
      allocate(sa(nfreq2, k))
      sa = dble(abs(yk)**2)
      evp = 1.0D0 - ev2
      tol = 0.03D0
      maxadit = 100
      call mw2wta(sa, dk, nfreq2, k, spec, ev2, evp, dofs, dofav, var &
        , dt, tol, maxadit, mxiter, aviter)
    end if

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

  subroutine weighted_eigencoef(d, ndata, dt, nw, k, yk, spec, nFFT, id)
  ! calculates and returns the weighted eigencoefficients
  !
  ! d   - real*8(ndata) - time series (data)
  ! ndata - integer     - length of the time series data
  ! dt    - real*8      - sampling rate
  ! nw    - real*8      - time bandwidth parameter
  ! k     - integer     - number of tapers to use
  ! yk - complex*16(nfreq, k) - matrix to hold the weighted eigencoefficients
  ! nFFT  - integer     - number of frequency bins to use (zeropadding)
  ! id    - integer     - which data set to use (currently 1 or 2)
  !
    integer :: ndata, nFFT, k, id
    real*8 :: d(ndata), nw, dt, var, dofav, spec(:)
    real*8, allocatable :: dk(:, :), dofs(:)
    complex*16 :: yk(:, :)

    if (id .eq. 1) then
      allocate(dk(nfreq, k), dofs(nfreq))
    else if (id .eq. 2) then
      allocate(dk(nfreq2, k), dofs(nfreq2))
    end if

    call variance(d, ndata, var)
    call eigencoef(d, ndata, nw, k, yk, nFFT, dt, id)
    call adaptive_weights(yk, dk, k, spec, dofs, dofav, var, dt, id)

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

  subroutine findHidx(col, hIdx, ncol, idxSub, baseIdx)
    ! (col, hIdx, nCol, nhIdx, idxSub, nSub) <- old arguments
    ! integer :: col(nCol), hIdx(nhIdx), idxAll(nCol), idxSub(nsub) &
    !  , ncol, nhIdx, nSub, i, j <- old integer declaration
  ! col - contains the indicator column of MSC matrix
  ! hIdx - contains the indices of all the unique frequencies used in the
  ! regression
  ! idxSub - contains the index of the column of H in which to store coefs
    integer :: col(:), hIdx(:), idxAll(ncol), idxSub(:), baseIdx &
      , ncol, i, j

    idxAll = (/ (i, i = 1, nCol) /)
    idxSub = pack(idxAll, col > 0 )
    j = 1

    do i = 1, size(idxSub)
      do j = j, size(hIdx)
        if (idxSub(i) == hIdx(j)) then
          idxSub(i) = j + baseIdx
          exit
        end if
      end do
    end do
  end subroutine findHidx


  subroutine msc_indicator(msc, nrow, ncol, ind, level, nOff)
  ! determines if the mcs's are local maxes (by column, i.e., central freq)
  ! and if the value exceeds the level provided (calculated based on the
  ! normal distribution ... stuff)
    real*8 :: msc(:, :), cohLmax(nrow, ncol), level
    integer :: i, j, ind(:, :), nrow, ncol, nOff, cmax(ncol)

    ind(:, :) = 0
    cohLmax(:, :) = -1

  ! only keep the local maxes above the cutoff.
    do i = 2, nrow-1
      do j = 1, ncol
        if (msc(i,j) >= level .and. msc(i, j) > msc(i-1, j) .and. &
          msc(i, j) > msc(i+1, j)) then
          cohLmax(i, j) = msc(i, j)
        end if
      end do
    end do

    do i = 1, nOff
      cmax = maxloc(msc, dim = 1, mask = (ind .eq. 0 .and. cohLmax > 0))
      do j = 1, ncol
        if (cmax(j) .eq. 0) cycle
        ind(cmax(j), j) = i
      end do
    end do
  end subroutine msc_indicator

  subroutine numberOffsets(col_order, total_offsets, totFreqByCol &
    , nrow, ncol, npred)
    integer :: total_offsets, totFreqByCol(ncol), wrktot(npred) &
      , wrk3(nrow, ncol), i, nrow, ncol, npred &
      , col_order(nrow, ncol, npred), wrk(nrow, npred)

      ! determine the number of predictors, including offsets (in total needed)
        ! number of total unique frequencies
        totFreqByCol = 0

        do i = 1, npred
          ! vector indicating which rows (offsets) have values to use
          wrk(:, i) = any(col_order(:, :, i) > -1, 2)
          ! the total number of unique offsets (per predictor)
          wrktot(i) = sum(wrk(:, i))
          total_offsets = total_offsets + wrktot(i)

        ! gets info for each design matrix
          ! indicator matrix - contains only 1's or 0's
          wrk3 = col_order(:, :, i) > -1
          ! total frequencies needed by center frequency - vector of length n_col
          totFreqByCol = totFreqByCol + sum(wrk3, 1)
        end do
  end subroutine

  subroutine setFreqIdx(freqIdx, fRatio)
    integer :: fRatio, freqIdx(nfreq2), i
    freqIdx = (/ (i, i = 1, nfreq2, fRatio) /)
  end subroutine setFreqIdx

  subroutine msc2norm(msc, nrow, ncol, dof)
    integer :: nrow, ncol, i, j, dof
    real*8 :: msc(nrow, ncol), tran

    do j = 1, ncol
      do i = 1, nrow
        tran = 1 - (1 - msc(i,j) )**(dof - 1)
        msc(i, j) = dsqrt(2.0D0) * erfinv(2*tran - 1)
      end do
    end do
  end subroutine msc2norm

  function cabs2(x, nrow, ncol)
    integer :: nrow, ncol
    real*8 :: cabs2(nrow, ncol)
    complex*16 :: x(nrow, ncol)

    cabs2 = realpart(x)**2 + imagpart(x)**2
  end function cabs2

  function erfinv(x)
    real*8 :: erfinv, x
    x = (x + 1.0D0) / 2.0D0
    erfinv = dinvnorm(x) / dsqrt(2.0D0)
  end function erfinv

  subroutine testem(x)
    integer :: x

    x = 5
  end subroutine testem

! taken from WayBackMachine: https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/
! and:
! https://stackedboxes.org/2017/05/01/acklams-normal-quantile-function/
  real*8 function dinvnorm(p)
  ! returns the inverse CDF at probability p
  ! i.e., the same as qnorm(p, 0, 1) in R
    ! real*8 :: dinvnorm
    real*8 :: p,p_low,p_high
    real*8 :: a1,a2,a3,a4,a5,a6
    real*8 :: b1,b2,b3,b4,b5
    real*8 :: c1,c2,c3,c4,c5,c6
    real*8 :: d1,d2,d3,d4
    real*8 :: z,q,r
    a1=-39.6968302866538
    a2=220.946098424521
    a3=-275.928510446969
    a4=138.357751867269
    a5=-30.6647980661472
    a6=2.50662827745924
    b1=-54.4760987982241
    b2=161.585836858041
    b3=-155.698979859887
    b4=66.8013118877197
    b5=-13.2806815528857
    c1=-0.00778489400243029
    c2=-0.322396458041136
    c3=-2.40075827716184
    c4=-2.54973253934373
    c5=4.37466414146497
    c6=2.93816398269878
    d1=0.00778469570904146
    d2=0.32246712907004
    d3=2.445134137143
    d4=3.75440866190742
    p_low=0.02425
    p_high=1-p_low
    if(p.lt.p_low) goto 201
    if(p.ge.p_low) goto 301
  201   q=dsqrt(-2*dlog(p))
    z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1)
    goto 204
  301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
    if(p.gt.p_high) goto 302
  202   q=p-0.5
    r=q*q
    z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
    goto 204
  302   if((p.gt.p_high).and.(p.lt.1)) goto 203
  203   q=dsqrt(-2*dlog(1-p))
    z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1)
  204   dinvnorm=z
    return
  end function dinvnorm

  integer function block_increment(block_size, overlap)
  ! determines index change for blocking based on overlap and block size.
    integer :: block_size
    real*8 :: overlap

    block_increment = floor(block_size * (1.0D0 - overlap))
    return
  end function block_increment

  integer function calculate_nblocks(ndata, block_size, block_incr)
  ! how many blocks given ndata data points, with block_size and
  ! increment.
    integer :: ndata, block_size, block_incr, j

    calculate_nblocks = size( (/ (j, j = 1 &
      , ndata - block_size + 1, block_incr) /) )
    return
  end function calculate_nblocks

  subroutine calc_tf_wt_eigen(block_incr, block_incr2 &
    , block_size, block_size2, nblocks, d1, d2, npred &
    , dt, dt2, nw, nw2, k, nFFT, nFFT2, yk1, yk2)
  ! does all the weighted eigencoefficient calculations and declutters
  ! tf() - also makes this more testable.
  ! must call fft_setup for series 1 and 2 before calling this subroutine

    integer :: block_incr, block_incr2, block_size, block_size2, nblocks &
      , k, nFFT, nFFT2, i, j, dstart_idx, dend_idx &
      , dstart_idx2, dend_idx2, npred
    real*8 :: d1(:), d2(:, :), dt, dt2, nw, nw2, s1(nfreq), s2(nfreq2)
    complex*16 :: yk1(:, :, :), yk2(:, :, :, :)

    do i = 0, nblocks - 1
      dstart_idx = i*block_incr + 1
      dend_idx = dstart_idx + block_size - 1
      dstart_idx2 = i*block_incr2 + 1
      dend_idx2 = dstart_idx2 + block_size2 - 1
      ! store the eigencoefficients by block for the response
      call weighted_eigencoef(d1(dstart_idx:dend_idx), block_size, dt &
       , nw, k, yk1(:, :, i+1), s1, nFFT, 1)
    ! store the eigencoefficients by block and predictor for the inputs
      do j = 1, npred
        call weighted_eigencoef(d2(dstart_idx2:dend_idx2, j) &
          , block_size2, dt2, nw2, k, yk2(:, :, j, i+1), s2, nFFT2, 2)
      end do
    end do
  end subroutine calc_tf_wt_eigen

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Testing section !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mtmtstcalctfwteigen(block_incr, block_incr2 &
    , block_size, block_size2, nblocks, d1, d2, dt, dt2, nw, nw2, k &
    , nFFT, nFFT2, yk1, yk2, ndata, ndata2, npred)

    integer :: block_incr, block_incr2, block_size, block_size2, nblocks &
      , k, nFFT, nFFT2, i, j, dstart_idx, dend_idx &
      , dstart_idx2, dend_idx2, ndata, ndata2, npred
    real*8 :: d1(ndata), d2(ndata2, npred), dt, dt2, nw, nw2, s1(nfreq), s2(nfreq2)
    complex*16 :: yk1(nfreq, k, nblocks), yk2(nfreq2, k, npred, nblocks)

    call calc_tf_wt_eigen(block_incr, block_incr2 &
      , block_size, block_size2, nblocks, d1, d2, npred &
      , dt, dt2, nw, nw2, k, nFFT, nFFT2, yk1, yk2)
  end subroutine mtmtstcalctfwteigen


  subroutine tstfindhidx(col, hIdx, ncol, nhIdx, idxSub, nSub, baseIdx)
    implicit none

    integer :: ncol, nhIdx, nSub, baseIdx, col(ncol), hIdx(nhIdx), idxSub(nSub)

    call findHidx(col, hIdx, ncol, idxSub, baseIdx)

  end subroutine tstfindhidx
end module mtm_mod
