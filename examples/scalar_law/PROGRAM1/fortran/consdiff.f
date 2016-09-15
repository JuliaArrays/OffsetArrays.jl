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
      subroutine consdiff(fic,lac,fif,laf,fix,lax,ifirst,ilast,
     &  x,flux,
     &  conservd)
      integer fic,lac,fif,laf,fix,lax,ifirst,ilast
      double precision
     &  x(fix:lax+1),
     &  flux(fif:laf+1)
      double precision
     &  conservd(fic:lac)
      integer ic
c     ******************************************************************
c     update conservd to new time
c     ******************************************************************
      do ic=ifirst,ilast
        conservd(ic) = conservd(ic)
     &               - (flux(ic+1)-flux(ic)) / (x(ic+1)-x(ic))
      enddo
      return
      end

      function stabledt(fc,lc,fm,lm,ifirst,ilast,
     &  conserved,lambda_cell,x)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "./const.i"
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision stabledt
      integer fc,lc,fm,lm,ifirst,ilast
      double precision conserved(fc:lc),lambda_cell(fc:lc),x(fm:lm+1)
      integer ic
      double precision abs_lambda,dx
c     ******************************************************************
c     compute stable timestep
c     ******************************************************************
      stabledt=huge
      do ic=ifirst,ilast
        abs_lambda=abs(lambda_cell(ic))
        dx = x(ic+1)-x(ic)
        if (abs_lambda.gt.roundoff*dx) then
          stabledt=min(stabledt,dx/abs_lambda)
        endif
      enddo
c     if (stabledt.le.0.) call abort()
      return
      end 
