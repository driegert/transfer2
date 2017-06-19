subroutine callblockcoh(d1, d2, ndata, block_size, overlap, dt &
  , nw, k, nFFT, coh, cohnrow, cohncol, freq, offsets, freq_range, max_freq_offset, calc_type, is_forward)
  use mtm_mod
  implicit none
  integer :: ndata, block_size, k, nFFT, calc_type, is_forward, cohnrow, cohncol
  real*8 :: d1(ndata), d2(ndata), overlap, dt, nw &
    , freq_range(2), max_freq_offset, freq(nFFT/2+1), offsets(cohnrow)
  complex*16 :: coh(cohnrow, cohncol)

  call block_coh(d1, d2, ndata, block_size, overlap, dt &
    , nw, k, nFFT, coh, freq, offsets, freq_range, max_freq_offset &
    , calc_type, is_forward)
end subroutine callblockcoh
