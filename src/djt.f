C$$$    The multitaper R package
C$$$    Multitaper and spectral analysis package for R
C$$$    Copyright (C) 2011 Karim J. Rahim David J. Thomson
C$$$
C$$$    This file is part of the multitaper package for R.
C$$$
C$$$    The multitaper package is free software: you can redistribute it and
C$$$    or modify
C$$$    it under the terms of the GNU General Public License as published by
C$$$    the Free Software Foundation, either version 2 of the License, or
C$$$    any later version.
C$$$
C$$$    The multitaper package is distributed in the hope that it will be
C$$$    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C$$$    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C$$$    GNU General Public License for more details.
C$$$
C$$$    You should have received a copy of the GNU General Public License
C$$$    along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
C$$$
C$$$    If you wish to report bugs please contact the author.
C$$$    karim.rahim@gmail.com
C$$$

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc This files contains modified djt multitaper files originally
cc from libraries written at Bell Labs by David Thomson.

c**************************************************************************
c     mw2wta multitaper using weights

      subroutine mw2wta(sa,wt,nfreq,nord,s,ev,evp
     1     ,dofs,dofav,var,dt,tol, maxadit
     1     , mxiter, aviter)
! sa    - real*8(nfreq, nord) - eigenspectra
! wt    - real*8(nfreq, nord) - adativeweights ... I think?
! nfreq - integer - nFFT/2 + 1 (half of the frequency bins)
! nord  - integer - number of tapers
! s     - real*8(nfreq) - adaptively weighted spectrum (+'ve freq)
! ev    - real*8(nord) - eigenvalues
! evp   - real*8(nord) - 1 - ev
! dofs  - real*8(nfreq) - dofs (approximate?)
! dofav - real*8(nfreq) - average dofs ? (I guess..)
! var   - real*8 - time series variance
! dt    - real*8 - sampling rate (?)
! tol   - real*8 - when to stop iterations
! maxadit - integer - maximum adaptive iterations (100 in R)
! mxiter  - integer - I don't know what this is
! aviter  - integer - average number of iterations?

      implicit none
      integer nfreq, nord, maxadit, mxiter,n, niter, k
      double precision  sa(nfreq,nord),wt(nfreq,nord),dofs(nfreq)
     1     ,s(nfreq),ev(nord),evp(nord),var,dt,tol, aviter, avewt
     1     ,wtmin, dofmin, sewn, sbar, dk2, sum, cwt, dofav, wmin
     1     ,dk2l


c     Generate Weights
       mxiter = 0
       aviter = 0.d0
       avewt = 0.d0
       wtmin = 1.d0
       dofmin = dble(2*nord)
       cwt = 0.d0
       wmin = 0.d0

c     Equivalent white noise level for bias calculation
       sewn = var*dt
      do 265 n=1,nfreq
c        start at estimate based on two best eigenvalues
         sbar = ( sa(n,1) + sa(n,2) )/2.d0
         dk2 = 1.d0
c     iterate
         do 262 niter=1, maxadit
            sum = 0.d0
            cwt = 0.d0
            wmin = 1.d0
            dk2l = dk2
            do 250 k = 1,nord
               dk2 =
     1              ( ev(k)*sbar/( ev(k)*sbar + evp(k)*sewn ) )**2
               wt(n,k) = dk2
               sum = sum + sa(n,k)*dk2
               wmin = dmin1(wmin,dk2)
               cwt = cwt + dk2
 250        continue
            sbar = sum/cwt
            if(dabs((dk2-dk2l)/(dk2+dk2l)).le.tol) exit
 262  continue
         mxiter = max0(mxiter,niter)
         aviter = aviter + niter
         avewt = avewt + cwt
         wtmin = dmin1(wtmin,wmin)
         dofs(n) = 2.d0*cwt
         dofmin = dmin1(dofmin,dofs(n))
         s(n) = sbar
         aviter = aviter/dble(nfreq)
 265  continue
      dofav = 2.d0*avewt/dble(nfreq)

      end subroutine
c*****end mw2wta

cc**********************************************************************
c     multiwindow jacknifed.

c     Multi-Window Weighting, Jackknifed
      subroutine mw2jkw(sa,wt,nfreq,nord,s,ev,evp
     1     ,dofs,dofav,var,dt,tol, sjk,varjk,bcjk,wjk,cwjk,vwj
     1     ,maxadit, mxiter)

      implicit none
      integer nfreq, nord, mxiter, n1, n2,j, n, niter, ks, maxadit
     1     ,k
      double precision sa(nfreq,nord),wt(nfreq,nord),dofs(nfreq)
     1     , s(nfreq)
     1     ,ev(nord),evp(nord),sjk(nord+2),varjk(nfreq),bcjk(nfreq)
     2     ,wjk(nord,nord+2),cwjk(nord+2),vwj(nord)
     1     ,total, avewt, wtmin, dofmin, bcor
     1     ,fnord, vnrm, sewn,var,dt, dofav, varcwt, sbar,sum
     1     ,wmin, slast, tol

c     Generate Weights

      mxiter = 0
      total = 0.d0
      avewt = 0.d0
      wtmin = 1.d0
      niter = 0
      sbar = 0.d0
      wmin =0.d0
      dofmin = dble(2*nord)
      bcor = dble(nord-1)
      fnord = nord
      vnrm = dble(nord-1)/fnord
      n1 = nord + 1
      n2 = nord + 2
c     Equivalent white noise level for bias calculation
      sewn = var*dt
c
      do 365 n=1,nfreq
c     iterate
         do 433 ks = 1, nord+1
c     start at estimate based on two best eigenvalues
            sbar = ( sa(n,1) + sa(n,2) )/2.
            do 362 niter=1, maxadit
               sum = 0.
               cwjk(ks) = 0.d0
               wmin = 1.d0
               slast = sbar
               do 350 k = 1,nord
                  if(k.eq.ks) go to 350
                  wjk(k,ks) = ( ev(k)*sbar/
     1                 ( ev(k)*sbar + evp(k)*sewn ) )**2
                  sum = sum + sa(n,k)*wjk(k,ks)
                  wmin = dmin1(wmin,wjk(k,ks))
                  cwjk(ks) = cwjk(ks) + wjk(k,ks)
 350           continue
               sbar = sum/cwjk(ks)
               sjk(ks) = dlog(sbar)
               if(dabs((sbar-slast)/(sbar+slast)).le.tol) exit
 362        continue
 433     continue

c     Jackknife mean, variance of Log S
         sjk(n2) = 0.d0
         cwjk(n2) = 0.d0
         do 490 k = 1, nord
            wjk(k,n2) = 0.d0
 490     continue
         do 500 k = 1, nord
            cwjk(n2) = cwjk(n2) + cwjk(k)
            sjk(n2) = sjk(n2) + sjk(k)
            do 510 j = 1, nord
               wjk(j,n2) = wjk(j,n2) + wjk(j,k)
 510        continue
 500     continue
         sjk(n2) = sjk(n2)/fnord
         cwjk(n2) = cwjk(n2)/fnord
         do 610 j = 1, nord
            vwj(j) = 0.d0
            wjk(j,n2) = wjk(j,n2)/fnord
            wt(n,j) = wjk(j,n2)
 610     continue

c     Jackknife Bias Estimate (Log S )
         bcjk(n) = bcor*(sjk(n2) - sjk(n1))

c     Variance Estimate
         varjk(n) = 0.d0
         varcwt = 0.d0
         do 550 k = 1, nord
            varjk(n) = varjk(n) + (sjk(k)-sjk(n2))**2
            varcwt = varcwt + (cwjk(k)-cwjk(n2))**2
            do 560 j = 1, nord
               vwj(j) = vwj(j) + (wjk(j,k)-wjk(j,n2))**2
 560        continue
 550     continue

         varjk(n) = varjk(n)*vnrm
         mxiter = max0(mxiter,niter)
         total = total + niter
         avewt = avewt + cwjk(n1)
         wtmin = dmin1(wtmin,wmin)
         dofs(n) = 2.d0*cwjk(n1)
         dofmin = dmin1(dofmin,dofs(n))
         s(n) = sbar
 365  continue
      dofav = 2.d0*avewt/float(nfreq)

      end subroutine
c ****** end mw2jkw


c     Multi-Window Average Estimation
      subroutine mweave(x,dw,swz,ndata,nord,ssqswz,cntr,dt
     1     ,spz, varc)
      implicit none
      integer ndata, nord, n, k, nnx
      double precision x(ndata),dw(ndata,nord),swz(nord)
     1     ,sm(nord),sum,spz, zero8, dt
     1     ,ssqswz, cntr, varc
      data  zero8/0.d+00/

c     no need for a max of 9
      nnx = nord

      call setdp(nnx,zero8,sm)
      do 100 k = 1, nnx
         do 110 n = 1, ndata
            sm(k) = sm(k) + dw(n,k)*x(n)
 110     continue
 100  continue

      sum = zero8
      spz = zero8
      do 300 k = 1, nnx, 2
         sum = sum + swz(k)*sm(k)
 300  continue
      sum = sum/ssqswz
      do 500 k = 1, nnx
         spz = spz + (sm(k) - sum*swz(k))**2
 500  continue
      spz = spz/dble(nnx)
      varc = spz/(dt*dble(ndata))
      cntr = sum
      end subroutine
c ******** end    mweave

c     Set Real*8 array
      subroutine setdp(npts,val,x)
      implicit none
      integer npts, n
      double precision val, x(npts)
      do 100 n = 1, npts
         x(n) = val
 100  continue

      end subroutine
c ******* end setdp

c ************************** helper functions used in coherence calculation
c djt/ts/ adstoa.f  Add Scalar to Array
      subroutine adstoa(x,y,ndata,xinc)
      implicit none
      integer n, ndata
      double precision  x(ndata), y(ndata), xinc
c     djt/ts/tsu1 -2- add scalar to array
      do 3400 n = 1,ndata
         y(n) = x(n)+xinc
 3400 continue

      end subroutine
c ********** end adstoa

c djt/ts/sphsed.f   Basic Phase Unwrapping Routine, Degrees
      subroutine sphsed(ph,nfreq)
      implicit none
      integer nfreq, n
      double precision ph(nfreq), q, pinc,d, t

      q=0.d0
      pinc=0.d0
      do 2100 n=1,nfreq
         t=ph(n)
         d=q-t
         q=t
         if(dabs(d).gt.180.d0) pinc=pinc+dsign(360.d0,d)
         ph(n)=t+pinc
 2100 continue

      end subroutine
c ****** end       sphsed


c*********************************************************************
