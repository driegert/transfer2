subroutine callblockcoh(d1, d2, ndata, ndata2, block_size, block_size2 &
  , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio, coh &
  , cohnrow, cohncol, freq, offsets, freq_range_idx &
  , max_freq_offset_idx, calc_type, conv_msc2norm, is_forward)

  use mtm_mod
  implicit none
  integer :: ndata, ndata2, block_size, block_size2, k, nFFT, nFFT2 &
    , fRatio, calc_type, is_forward, cohnrow, cohncol &
    , freq_range_idx(2), max_freq_offset_idx, conv_msc2norm
  real*8 :: d1(ndata), d2(ndata2), overlap, dt, dt2, nw, nw2 &
    , freq(nFFT/2+1), offsets(cohnrow)
  complex*16 :: coh(cohnrow, cohncol)

  call block_coh(d1, d2, ndata, ndata2, block_size, block_size2 &
    , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio, coh &
    , freq, offsets, freq_range_idx, max_freq_offset_idx, calc_type &
    , conv_msc2norm, is_forward)
end subroutine callblockcoh


subroutine calltf(d1, d2, ndata, ndata2, npred, block_size, block_size2 &
  , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio &
  , freq_range_idx, max_freq_offset_idx, H &
  , coh_nrow, coh_ncol, totFreqByCol, totFreqByColByPred, total_offsets &
  , col_order, hPredBreak, hIdx, nhIdx)
  use mtm_mod
  implicit none

  integer :: ndata, ndata2, npred, block_size, block_size2, k &
    , nFFT, nFFT2, fRatio, freq_range_idx(2), max_freq_offset_idx &
    , coh_nrow, coh_ncol, totFreqByCol(coh_ncol), total_offsets &
    , col_order(coh_nrow, coh_ncol, npred), hPredBreak(npred+1) &
    , hIdx(nhIdx), nhIdx, totFreqByColByPred(npred+1, coh_ncol)
  real*8 :: d1(ndata), d2(ndata2), overlap, dt, dt2, nw, nw2
  complex*16 :: H(coh_ncol, total_offsets)

  call tf(d1, d2, ndata, ndata2, npred, block_size, block_size2 &
    , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio &
    , freq_range_idx, max_freq_offset_idx, H &
    , coh_nrow, coh_ncol, totFreqByCol, totFreqByColByPred, total_offsets &
    , col_order, hPredBreak, hIdx, nhIdx)
end subroutine calltf


subroutine calltfzero(d1, d2, ndata, ndata2, npred, block_size, block_size2 &
  , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio &
  , freq_range_idx, H, n_row_H)
  use mtm_mod
  implicit none

  integer :: ndata, ndata2, npred, block_size, block_size2, k, nFFT, nFFT2 &
    , fRatio, freq_range_idx(2), n_row_H
  real*8 :: d1(ndata), d2(ndata, npred), overlap, dt, dt2, nw, nw2
  complex*16 :: H(n_row_H, npred)

  call tf_zero(d1, d2, ndata, ndata2, npred, block_size, block_size2 &
    , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, fRatio &
    , freq_range_idx, H, n_row_H)
end subroutine calltfzero

subroutine callMscIndicator(msc, nrow, ncol, ind, level, nOff)
  use mtm_mod
  implicit none

  integer :: nrow, ncol, nOff, ind(nrow, ncol)
  real*8 :: msc(nrow, ncol), level

  call msc_indicator(msc, nrow, ncol, ind, level, nOff)
end subroutine callMscIndicator
