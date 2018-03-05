subroutine tstWork(x)
  use mtm_mod
  integer :: x

  x = 54
end subroutine tstWork


! Different answers due to scaling by dt in different spots
subroutine tstTaperData(d1, d2, n1, n2, m1, m2, c1, c2, dt1, dt2, k, nw)
  use mtm_mod
  integer :: n1, n2, m1, m2, k
  real*8 :: d1(n1), d2(n2), dt1, dt2, nw
  complex*16 :: c1(m1, k), c2(m2, k)

  call fft_setup(m1, 1)
  call fft_setup(m2, 2)
  call dpss_setup(n1, nw, k, m1, 1)
  call dpss_setup(n2, nw, k, m2, 2)

  call taper_data(d1, c1, n1, k, m1, dt1, 1)
  call taper_data(d2, c2, n2, k, m2, dt2, 2)

  call dpss_cleanup(1)
  call dpss_cleanup(2)
  call fft_cleanup(1)
  call fft_cleanup(2)

end subroutine tstTaperData

subroutine tstEigenCoef(d1, d2, n1, n2, m1, m2, yk1, yk2, dt1, dt2, k, nw)
  use mtm_mod

  integer :: n1, n2, m1, m2, k
  real*8 :: d1(n1), d2(n2), dt1, dt2, nw
  complex*16 :: yk1( (m1/2)+1, k), yk2( (m2/2)+1, k)

  call fft_setup(m1, 1)
  call fft_setup(m2, 2)
  call dpss_setup(n1, nw, k, m1, 1)
  call dpss_setup(n2, nw, k, m2, 2)

  call eigenCoef(d1, n1, nw, k, yk1, m1, dt1, 1)
  call eigenCoef(d2, n2, nw, k, yk2, m2, dt2, 2)

  call dpss_cleanup(1)
  call dpss_cleanup(2)
  call fft_cleanup(1)
  call fft_cleanup(2)
end subroutine tstEigenCoef

subroutine tstWeightedEigenCoef(d1, d2, n1, n2, m1, m2, yk1, yk2, dt1, dt2, k, nw)
  use mtm_mod

  integer :: n1, n2, m1, m2, k
  real*8 :: d1(n1), d2(n2), dt1, dt2, nw, s1( (m1/2)+2 ), s2( (m2/2)+1 )
  complex*16 :: yk1( (m1/2)+1, k), yk2( (m2/2)+1, k)

  call fft_setup(m1, 1)
  call fft_setup(m2, 2)
  call dpss_setup(n1, nw, k, m1, 1)
  call dpss_setup(n2, nw, k, m2, 2)

  call weighted_eigencoef(d1, n1, dt1, nw, k, yk1, s1, m1, 1)
  call weighted_eigencoef(d2, n2, dt2, nw, k, yk2, s2, m2, 2)

  call dpss_cleanup(1)
  call dpss_cleanup(2)
  call fft_cleanup(1)
  call fft_cleanup(2)

end subroutine tstWeightedEigenCoef

! subroutine tstAdaptiveWeights(d1, d2, n1, n2, m1, m2, yk1, yk2 &
!   , dt1, dt2, k, nw, dk1, dk2, var1, var2)
!   ! this crashes for some reason... not *quite* sure why ...
!   use mtm_mod
!
!   integer :: n1, n2, m1, m2, k
!   real*8 :: d1(n1), d2(n2), dt1, dt2, nw, s1( (m1/2)+2 ), s2( (m2/2)+1 ) &
!     , dk1( (m1/2)+1, k), dk2( (m2/2)+1, k), dofs1((m1/2)+1), dofs2((m2/2)+1) &
!     , dofav1, dofav2, var1, var2
!   complex*16 :: yk1( (m1/2)+1, k), yk2( (m2/2)+1, k)
!
!   call fft_setup(m1, 1)
!   call fft_setup(m2, 2)
!   call dpss_setup(n1, nw, k, m1, 1)
!   call dpss_setup(n2, nw, k, m2, 2)
!
!   call eigenCoef(d1, n1, nw, k, yk1, m1, dt1, 1)
!   call eigenCoef(d2, n2, nw, k, yk2, m2, dt2, 2)
!
!   ! runs if these next two lines are commented out ... no clue why
!   call adaptive_weights(yk1, dk1, k, s1, dofs1, dofav1, var1, dt1, 1)
!   ! call adaptive_weights(yk2, dk2, k, s2, dofs2, dofav2, var2, dt2, 2)
!
!   call dpss_cleanup(1)
!   call dpss_cleanup(2)
!   call fft_cleanup(1)
!   call fft_cleanup(2)
!
! end subroutine tstAdaptiveWeights

subroutine tstEigenvals(n1, n2, k, nw, m1, m2, lambda1, lambda2)
  use mtm_mod
  implicit none
  integer :: n1, n2, k, m1, m2
  real*8 :: nw, lambda1(k), lambda2(k)


  call fft_setup(m1, 1)
  call fft_setup(m2, 2)
  call dpss_setup(n1, nw, k, m1, 1)
  call dpss_setup(n2, nw, k, m2, 2)
  call dpssToEigenvalues(n1, k, nw, m1, 1)
  call dpssToEigenvalues(n2, k, nw, m2, 2)
  lambda1 = ev
  lambda2 = ev2
  call dpss_cleanup(1)
  call fft_cleanup(1)
end subroutine tstEigenvals


subroutine tsthfind(col, hIdx, ncol, nhIdx, idxSub, nSub, baseIdx)
  use mtm_mod
  implicit none

  integer :: ncol, nhIdx, nSub, baseIdx, col(ncol), hIdx(nhIdx) &
    , idxSub(nSub)

  call tstfindhidx(col, hIdx, ncol, nhIdx, idxSub, nSub, baseIdx)
end subroutine tsthfind


subroutine tstblockincrement(block_size, overlap, block_incr)
  use mtm_mod
  implicit none

  integer :: block_size, block_incr
  real*8 :: overlap

  block_incr = block_increment(block_size, overlap)
end subroutine tstblockincrement


subroutine tstcalculatenblocks(ndata, block_size, block_incr, nblocks)
  use mtm_mod
  implicit none

  integer :: ndata, block_size, block_incr, nblocks

  nblocks = calculate_nblocks(ndata, block_size, block_incr)
end subroutine tstcalculatenblocks

subroutine tstcalctfwteigen(block_incr, block_incr2 &
  , block_size, block_size2, nblocks, d1, d2, dt, dt2, nw, nw2, k &
  , nFFT, nFFT2, yk1, yk2, ndata, ndata2, npred)
  use mtm_mod
  implicit none

  integer :: block_incr, block_incr2, block_size, block_size2, nblocks &
    , k, nFFT, nFFT2, i, j, dstart_idx, dend_idx &
    , dstart_idx2, dend_idx2, ndata, ndata2, npred, nrow, ncol
  real*8 :: d1(ndata), d2(ndata2, npred), dt, dt2, nw, nw2, s1(nfreq), s2(nfreq2)
  complex*16 :: yk1(nfreq, k, nblocks), yk2(nfreq2, k, npred, nblocks)

  call fft_setup(nFFT, 1)
  call fft_setup(nFFT2, 2)
  call dpss_setup(block_size, nw, k, nFFT, 1)
  call dpss_setup(block_size2, nw2, k, nFFT2, 2)

  call mtmtstcalctfwteigen(block_incr, block_incr2 &
    , block_size, block_size2, nblocks, d1, d2, dt, dt2, nw, nw2, k &
    , nFFT, nFFT2, yk1, yk2, ndata, ndata2, npred)

  call fft_cleanup(1)
  call fft_cleanup(2)
  call dpss_cleanup(1)
  call dpss_cleanup(2)
end subroutine tstcalctfwteigen

subroutine tsttfzero(d1, d2, ndata, ndata2, npred, block_size, block_size2 &
  , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio &
  , freq_range_idx, H, n_row_H)
  use mtm_mod
  implicit none

  integer :: ndata, ndata2, npred, block_size, block_size2, k &
    , nFFT, nFFT2, fRatio, freq_range_idx(2), n_row_H
  real*8 :: d1(ndata), d2(ndata2, npred), nw, nw2, overlap, dt, dt2
  complex*16 :: H(n_row_H, npred)

  call tf_zero(d1, d2, ndata, ndata2, npred, block_size, block_size2 &
    , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio &
    , freq_range_idx, H, n_row_H)
end subroutine tsttfzero

subroutine tstmatrixwrap(x)
  use mtm_mod
  integer :: x(5,3)
  call tstmatrix(x)
end subroutine tstmatrixwrap
