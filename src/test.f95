subroutine adder(x)
  real*8 :: x
  x = 2*x +9
end subroutine adder

subroutine test2(x)
  use mtm_mod
  implicit none
  integer :: x

  call testem(x)
end subroutine test2

subroutine testcd(cs, s1, s2, cnrow, cncol &
  , idx_start, idx_end, max_idx, nFFT)
  use mtm_mod
  implicit none

  integer :: cnrow, cncol, idx_start, idx_end, max_idx, nFFT
  complex*16 :: cs(cnrow, cncol)
  real*8 :: s1(nFFT/2+1), s2(nFFT/2+1)

  call fft_setup(nFFT)
  call coh_denom(cs, s1, s2, cnrow, cncol, idx_start, idx_end, max_idx)
end subroutine testcd

subroutine test3(nf)
  use mtm_mod
  implicit none
  integer :: nf

  call fft_setup(128)
  call testem(nf)
end subroutine

subroutine testdpss(n, nw, k, vec, val, nFFT)
  use mtm_mod
  implicit none
  integer :: n, k, nFFT
  real*8 :: nw, vec(n, k), val(k)

  call fft_setup(nFFT)
  call dpss_setup(n, nw, k, nFFT)

  vec = v
  val = ev

  call fft_cleanup
  call dpss_cleanup
end subroutine testdpss


! test the fourier transforms for the eigenvalues ... see what the hell is going on.
subroutine testratios(ndata, nw, k, ratios, vec, val, nFFT, cx)
  use mtm_mod
  implicit none
  integer :: j, ndata, k, nFFT
  real *8 :: nw, w, ratios(ndata), vec(ndata, k), val(k)
  complex*16 :: cx(nFFT, k)

  call fft_setup(nFFT)
  call dpss_setup(ndata, nw, k, nFFT)

  vec = v
  val = ev

  w = nw / dble(ndata)

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
    ev(j) = 2.0D0*ev(j) + ratios(1)*realpart(cx(1, j))
  end do

  call fft_cleanup
  call dpss_cleanup
end subroutine testratios

subroutine testcross(yk1, yk2, cs12, k, nFFT, idx_start, idx_end, idx_max)
  use mtm_mod
  implicit none

  integer :: k, nFFT, idx_start, idx_end, idx_max, forward
  complex*16 :: yk1(nFFT/2+1, k), yk2(nFFT/2+1, k), cs12(2*idx_max + 1, idx_end - idx_start + 1)

  forward = 1

  call fft_setup(nFFT)
  call cross_spec(yk1, yk2, cs12, k, idx_start, idx_end, idx_max, forward)
  call fft_cleanup
end subroutine testcross

subroutine testeigen(d, ndata, nw, k, yk, nFFT, dt)
  use mtm_mod
  implicit none
  integer :: ndata, k, nFFT
  real*8 :: d(ndata), nw, dt, s(nFFT/2+1)
  complex*16 :: yk(nFFT, k)

  call fft_setup(nFFT)
  call dpss_setup(ndata, nw, k, nFFT)
  ! call eigenCoef(d, ndata, nw, k, yk, nFFT, dt)
  call weighted_eigencoef(d, ndata, dt, nw, k, yk, s, nFFT)

  call fft_cleanup
  call dpss_cleanup
end subroutine testeigen

subroutine testweights(d, yk, dk, k, var, dt, ndata, nFFT, nw, eval)
  use mtm_mod
  implicit none
  integer :: ndata, nFFT, k
  real*8 :: d(ndata), nw, dt, var, dk(nFFT/2+1, k), spec(nFFT/2+1), dofs(nFFT/2+1), dofav, eval(k)
  complex*16 :: yk(nFFT/2+1, k)

  call fft_setup(nFFT)
  call dpss_setup(ndata, nw, k, nFFT)

  ev = eval
  ! call variance(d, ndata, var)
  ! call eigencoef(d, ndata, nw, k, yk, nFFT, dt)
  call adaptive_weights(yk, dk, k, spec, dofs, dofav, var, dt)
  call fft_cleanup
  call dpss_cleanup
end subroutine testweights

subroutine weightwrapper(sa, dk, nFFT, k, ev, var, dt)
  implicit none
  integer :: k, nFFT, maxadit, mxiter, nfreq
  real*8 :: sa(nFFT/2+1, k), dofs(nFFT/2+1)
  real*8 :: dofav, dk(nFFT/2+1, k), spec(nFFT/2+1), ev(k), evp(k), tol, aviter, var, dt

  nfreq = nFFT/2+1
  evp = 1.0D0 - ev
  tol = 0.03D0
  maxadit = 100
  call mw2wta(sa, dk, nfreq, k, spec, ev, evp, dofs, dofav, var, dt, tol, maxadit, mxiter, aviter)
  dk = dsqrt(dk)

end subroutine weightwrapper

subroutine testfft(d, ndata, nFFT, c)
  use mtm_mod
  implicit none

  integer :: ndata, nFFT
  real*8 :: d(ndata)
  complex*16 :: c(nFFT)

  c(1:ndata) = dcmplx(d(1:ndata), 0.0D0)
  call fft_setup(nFFT)
  call cfft1f(nFFT, inc, c, lenc, wsave, lensav, work, lenwrk, ier)
  call fft_cleanup
end subroutine testfft

! subroutine testpfft(d, ndata, nFFT, c)
!   implicit none
!   integer :: ndata, nFFT
!   real*8 :: d(ndata), ar(nFFT), ai(nFFT)
!   complex*16 :: c(nFFT)
!
!   ar = 0.0D0
!   ai = 0.0D0
!   ar(1:ndata) = d(1:ndata)
!   call dfftc(nFFT, ar, ai)
!
!   c = dcmplx(ar, ai)
! end subroutine testpfft

subroutine testcross2(d1, d2, ndata, cs12, nw, k, nFFT, idx_start, idx_end, idx_max, dt)
  use mtm_mod
  implicit none

  integer :: ndata, k, nFFT, idx_start, idx_end, idx_max, forward
  real*8 :: d1(ndata), d2(ndata), nw, s1(nFFT/2+1), s2(nFFT/2+1), dt
  complex*16 :: yk1(nFFT/2+1, k), yk2(nFFT/2+1, k), cs12(2*idx_max + 1, idx_end - idx_start + 1)

  forward = 1

  call fft_setup(nFFT)
  call dpss_setup(ndata, nw, k, nFFT)
  call weighted_eigencoef(d1, ndata, dt, nw, k, yk1, s1, nFFT)
  call weighted_eigencoef(d2, ndata, dt, nw, k, yk2, s2, nFFT)
  call cross_spec(yk1, yk2, cs12, k, idx_start, idx_end, idx_max, forward)
  call fft_cleanup
  call dpss_cleanup
end subroutine testcross2

subroutine svdRegTest(Y, X, m, n, beta, stdErr, svd_ev)
  use mtm_mod
  implicit none

  integer :: m, n
  real*8 :: svd_ev(n), stdErr(n)
  complex*16 :: Y(m), X(m, n), beta(n)

  call zSvdRegression(Y, X, m, n, beta, stdErr, svd_ev)
end subroutine svdRegTest
