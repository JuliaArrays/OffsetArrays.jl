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
      subroutine initsl(ncells,fc,lc,fx,lx,ifirst,ilast,
     &  conservd,x)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "linearad.i"
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer ncells,fc,lc,fx,lx,ifirst,ilast
      double precision 
     &  conservd(fc:lc),
     &  x(fx:lx+1)
      integer ic,ie,ijump
      double precision dx,frac
c     ******************************************************************
c     interior values only; others defined by boundary conditions
c     ******************************************************************
      dx=(x_right-x_left)/dble(ncells)
      do ie=ifirst,ilast+1
        x(ie)=x_left+dble(ie)*dx
      enddo
      ijump=max(ifirst-1,min(
     &        int(dble(ncells)*(jump-x_left)/(x_right-x_left)),ilast+1))
c     print *, "fc,lc = ",fc,lc
c     print *, "ifirst,ilast = ",ifirst,ilast
c     print *, "jump,x_left,x_right = ",jump,x_left,x_right
      do ic=ifirst,ijump-1
        conservd(ic)=statelft
      enddo
      frac=(jump-x_left-dble(ijump)*dx)/(x_right-x_left)
      conservd(ijump)=statelft*frac+statergt*(1.d0-frac)
      do ic=ijump+1,ilast
        conservd(ic)=statergt
      enddo
      return
      end 

      subroutine bcmesh(fim,lam,ncells,
     &  x)
      integer fim,lam,ncells
      double precision  
     &  x(fim:lam+1)
      integer ie,it
      double precision dx
c     ******************************************************************
c     right side 
c     ******************************************************************
      if (lam.ge.ncells) then
        it=ncells
        dx=x(it)-x(it-1)
        do ie=it+1,lam+1
          x(ie)=x(it)+dx*dble(ie-it)
        enddo
      endif
c     ******************************************************************
c     left side 
c     ******************************************************************
      if (fim.lt.0) then
        dx=x(1)-x(0)
        do ie=fim,-1
          x(ie)=x(0)+dx*dble(ie)
        enddo
      endif
      return
      end

      subroutine bccells(fic,lac,ncells,
     &  conservd)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "linearad.i"
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fic,lac,ncells
      double precision 
     &  conservd(fic:lac)
      integer ic
c     ******************************************************************
c     right side 
c     ******************************************************************
      if (lac.ge.ncells) then
c       print *, "treating right side"
        do ic=ncells,lac
c         outgoing wave:
          conservd(ic)=conservd(ncells-1)
        enddo
      endif
c     ******************************************************************
c     left side 
c     ******************************************************************
      if (fic.lt.0) then
c       print *, "treating left side"
        do ic=fic,-1
c         outgoing wave:
c         conservd(ic)=conservd(0)

c         specified value
          conservd(ic)=statelft
        enddo
      endif
      return
      end
