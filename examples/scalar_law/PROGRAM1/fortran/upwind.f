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
      subroutine method(dt,fc,lc,fm,lm,fs,ls,ifirst,ilast,
     &  conserved,
     &  flux)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "linearad.i"
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fc,lc,fm,lm,fs,ls,ifirst,ilast
      double precision dt
      double precision conserved(fc:lc)
      double precision flux(fs:ls+1)
      integer ie
      double precision vdt
c     ******************************************************************
c     multiply fluxes by dt (ie, compute time integral over cell side)
c     ******************************************************************
      vdt=velocity*dt
      do ie=ifirst,ilast+1
        flux(ie)=vdt*conserved(ie-1)
      enddo

      return
      end
