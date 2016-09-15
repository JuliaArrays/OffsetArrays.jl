c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      program linaddmain
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer ncells
!      parameter (ncells=100)
!      parameter (ncells=10000)
      parameter (ncells=100000)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "linearad.i"
      include "const.i"
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision stabledt
      external stabledt

      double precision
     &  u(-2:ncells+1),
     &  x(0:ncells),
     &  flux(0:ncells),
     &  dfdu(-2:ncells+1)
c    &  flux_cell(-2:ncells+1),
c    &  uside(0:ncells)
      integer nsteps
      double precision tmax,cfl

      integer i,fc,lc,fm,lm,fs,ls,ifirst,ilast
      integer istep
      double precision dt,t
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      roundoff=1.d-14
      small=1.d-20
      huge=1.d300
      undefind=1.d300
      lnrndoff=14.d0
      lnsmall=20.d0

c     problem-specific parameters:
      jump=0
      x_left=-0.2d0
      x_right=1.d0
      statelft=2.d0
      statergt=0.d0
      velocity=1.d0

      nsteps=10000
      tmax=0.8d0
      cfl=0.9d0

c     array bounds for subroutine calls:
      fc=-2
      lc=ncells+1
      fm=0
      lm=ncells-1
      fs=0
      ls=ncells-1
      ifirst=0
      ilast=ncells-1

c     initialization:
      call initsl(ncells,fc,lc,fm,lm,ifirst,ilast, u,x)
c     do i=0,ncells
c       print *,"x(",i,") = ",x(i)
c     enddo
c     do i=-2,ncells+1
c       print *,"u(",i,") = ",u(i)
c     enddo
      call bcmesh(fm,lm,ncells, x)
      call bccells(fc,lc,ncells, u)
      call fluxderv(fc,lc,fc,lc, u, dfdu)

      istep=0
      t=0.d0
      do while (istep.lt.nsteps .and. t.lt.tmax)
        call bccells(fc,lc,ncells, u)
        call fluxderv(fc,lc,fc,lc, u, dfdu)
        dt=cfl*stabledt(fc,lc,fm,lm,ifirst,ilast, u, dfdu,x)
        call method(dt,fc,lc,fm,lm,fs,ls,ifirst,ilast, u, flux)
        call consdiff(fc,lc,fs,ls,fm,lm,ifirst,ilast, x,flux, u)
        t=t+dt
        istep=istep+1
      enddo

c     write results, plot later
      do i=0,ncells-1
!        print *,(x(i)+x(i+1))*0.5d0,u(i)
      enddo
     
      stop
      end
