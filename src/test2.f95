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
