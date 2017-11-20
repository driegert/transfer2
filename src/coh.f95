subroutine callblockcoh(d1, d2, ndata, ndata2, block_size, block_size2 &
  , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, coh, cohnrow, cohncol &
  , freq, offsets, freq_range_idx &
  , max_freq_offset_idx, calc_type, conv_msc2norm, is_forward)

  use mtm_mod
  implicit none
  integer :: ndata, ndata2, block_size, block_size2, k, nFFT, nFFT2 &
    , calc_type, is_forward, cohnrow, cohncol, freq_range_idx(2) &
    , max_freq_offset_idx, conv_msc2norm
  real*8 :: d1(ndata), d2(ndata2), overlap, dt, dt2, nw, nw2 &
    , freq(nFFT/2+1), offsets(cohnrow)
  complex*16 :: coh(cohnrow, cohncol)

  call block_coh(d1, d2, ndata, ndata2, block_size, block_size2 &
    , overlap, dt, dt2, nw, nw2, k, nFFT, nFFT2, coh, freq, offsets &
    , freq_range_idx, max_freq_offset_idx, calc_type &
    , conv_msc2norm, is_forward)
end subroutine callblockcoh


subroutine calltf()
  use mtm_mod
  implicit none

end subroutine calltf
