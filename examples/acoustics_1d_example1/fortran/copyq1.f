c     
c     
c=========================================================
      subroutine copyq1(meqn,mbc,mx,q1,q2)
c=========================================================
c     
c     # copy the contents of q1 into q2
c     
      implicit double precision (a-h,o-z)
      dimension q1(meqn,1-mbc:mx+mbc)
      dimension q2(meqn,1-mbc:mx+mbc)
c     
      do i = 1-mbc, mx+mbc
         do m=1,meqn
            q2(m,i) = q1(m,i)
         end do
      end do
      return
      end

