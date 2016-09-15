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
      subroutine fluxderv(fc,lc,ifirst,ilast,
     &  conservd,
     &  dfdu    )
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "const.i"
      include "linearad.i"
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fc,lc,ifirst,ilast
      double precision conservd(fc:lc)
c
      double precision dfdu(fc:lc)
c
      integer ic
c     ******************************************************************
c     compute derivative of flux with respect to conserved variable
c     ******************************************************************
      do ic=ifirst,ilast
        dfdu(ic)=velocity
      enddo

      return
      end 

      subroutine riemann(fs,ls,ifirst,ilast,
     &  left,right,
     &  flux)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "linearad.i"
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fs,ls,ifirst,ilast
c
      double precision  
     &  left(fs:ls+1),
     &  right(fs:ls+1)
c
      double precision flux(fs:ls+1)
      integer ie
c     ******************************************************************
c     the following riemann solver is hard-wired for linear advection
c     with velocity > 0
c     ******************************************************************
      do ie=ifirst,ilast+1
        flux(ie)=velocity*left(ie)
      enddo
      return
      end 
