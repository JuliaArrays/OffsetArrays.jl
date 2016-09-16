c     
c     
c=========================================================
      subroutine out1(meqn,mbc,mx,xlower,dx,q,t,iframe,aux,maux,outaux)
c=========================================================
c     
c     # Output the results for a general system of conservation laws
c     # in 1 dimension
c     
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw1.m
c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c     
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>
c     
c     
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      character*10 fname1, fname2, fname3
      logical outaux

c     
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw1.m
c     # The same format is used by the amrclaw package.  
c     # Here it's adapted to output just the single grid.
c     
c     # first create the file name and open file
c     
      fname1 = 'fort.qxxxx'
      fname2 = 'fort.txxxx'
      fname3 = 'fort.axxxx'
      nstp = iframe
      do 55 ipos = 10, 7, -1
         idigit = mod(nstp,10)
         fname1(ipos:ipos) = char(ichar('0') + idigit)
         fname2(ipos:ipos) = char(ichar('0') + idigit)
         fname3(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
 55   continue

      open(unit=50,file=fname1,status='unknown',form='formatted')
      open(unit=60,file=fname2,status='unknown',form='formatted')

c     
c     # the following parameters are used in amrclaw where there are
c     # multiple grids.  Here they are all set to 1:
      ngrids = 1
      mptr = 1
      level = 1

      write(50,1001) mptr,level,mx
 1001 format(i5,'                 grid_number',/,
     &     i5,'                 AMR_level',/,
     &     i5,'                 mx')

      write(50,1002) xlower,dx
 1002 format(e18.8,'    xlow', /,
     &     e18.8,'    dx', /)
c     
      do 10 i=1,mx
         do m=1,meqn
c     # exponents with more than 2 digits cause problems reading
c     # into matlab... reset tiny values to zero:
            if (dabs(q(m,i)) .lt. 1d-99) q(m,i) = 0.d0
         enddo
c     
         write(50,1005) (q(m,i), m=1,meqn)
 1005    format(4e16.8)
c     
 10   continue
      write(50,*) ' '
      write(50,*) ' '


      if (outaux) then 
c     # also output the aux arrays:
         open(unit=70,file=fname3,status='unknown',form='formatted')
         write(70,1001) mptr,level,mx
         write(70,1002) xlower,dx
         do 110 i=1,mx
            do m=1,maux
c     # exponents with more than 2 digits cause problems reading
c     # into matlab... reset tiny values to zero:
               if (dabs(aux(m,i)) .lt. 1d-99) aux(m,i) = 0.d0
            enddo
c     
            write(70,1005) (aux(m,i), m=1,maux)
c     
 110     continue
         write(70,*) ' '
         close(unit=70)
      endif

      write(60,1000) t,meqn,ngrids,maux,1

 1000 format(e26.16,'    time', /,
     &     i5,'                 meqn'/,
     &     i5,'                 ngrids'/,
     &     i5,'                 maux'/,
     &     i5,'                 ndim'/,/)
c     

      close(unit=50)
      close(unit=60)

      return
      end

