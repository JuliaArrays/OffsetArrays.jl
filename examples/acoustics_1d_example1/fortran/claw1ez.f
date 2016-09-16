c
c
c
c     ===============================================================
      subroutine claw1ez    ! No arguments
c     ===============================================================
c
c     An easy-to-use clawpack driver routine for simple applications
c     Documentation is available at
c                 http://www.amath.washington.edu/~claw/doc.html
c
c     Authors: Randall J. LeVeque, Grady I. Lemoine
c
      implicit double precision (a-h,o-z)
      external bc1,rp1,src1,b4step1

      double precision, dimension(:,:), allocatable :: q, aux
      double precision, dimension(:), allocatable :: work, tout
      integer, dimension(:), allocatable :: mthlim, iout_q, iout_aux

c
      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(2)
      integer :: allocate_status, outstyle
      logical :: outaux_init_only, use_fwaves, output_t0
      logical :: outaux_always, dt_variable
      character*12 fname
c
      open(10,file='fort.info',status='unknown',form='formatted')

c
c     # Read the input in standard form from claw.data:
c     # For a description of input parameters see the documentation at
c                 http://www.amath.washington.edu/~claw

c     ## New in 4.4:   Input file name changed from claw1ez.data
c     ## Open file and skip over leading lines with # comments:
      fname = 'claw.data'
      call opendatafile(55,fname)

c
c     # Read the input in standard form from claw.data:
c

      read(55,*) ndim

      read(55,*) xlower
      read(55,*) xupper
      read(55,*) mx
      read(55,*) meqn
      read(55,*) mwaves
      read(55,*) maux
      read(55,*) t0

      read(55,*) outstyle
      if (outstyle == 1) then
         read(55,*) nout
         read(55,*) tfinal
         read(55,*) output_t0    ! Not currently used
         nstepout = 1
      else if (outstyle == 2) then
         read(55,*) nout
         allocate(tout(nout), stat=allocate_status)
         if (allocate_status .ne. 0) then
            print *, '*** Error allocating tout array; exiting claw1ez'
            go to 900
         end if
         read(55,*) (tout(i), i=1,nout)
         nstepout = 1
      else if (outstyle == 3) then
         read(55,*) nstepout
         read(55,*) nstop
         read(55,*) output_t0
         nout = nstop
      else
         print *, '*** Unrecognized output style ', outstyle
         print *, '*** Exiting claw1ez'
         go to 900
      end if

      read(55,*) output_format    ! Not used yet
      ! These iout variables are not currently used, but hang onto them
      ! anyway in case somebody wants to use them at a future date.  The
      ! same goes for outaux_init_only.
      allocate(iout_q(meqn), stat=allocate_status)
      if (allocate_status .ne. 0) then
         print *, '*** Error allocating iout_q array; exiting claw1ez'
         go to 900    ! Exception handling, old school style
      end if
      read(55,*) (iout_q(i), i = 1, meqn)
      if (maux > 0) then
         allocate(iout_aux(maux), stat=allocate_status)
         if (allocate_status .ne. 0) then
            print *, '*** Error allocating iout_aux array;',
     &               ' exiting claw1ez'
            go to 900
         end if
         read(55,*) (iout_aux(i), i = 1, maux)
         read(55,*) outaux_init_only
         ! Not implementing selective output of aux fields yet
         if (any(iout_aux .ne. 0)) then
            outaux_always = .not. outaux_init_only
         else
            outaux_always = .false.
            outaux_init_only = .false.
         end if
      else
         outaux_always = .false.
         outaux_init_only = .false.    ! Just to initialize
      end if

      read(55,*) dtv(1)     ! Initial dt
      read(55,*) dtv(2)     ! Max dt
      read(55,*) cflv(1)    ! Max CFL number
      read(55,*) cflv(2)    ! Desired CFL number
      read(55,*) nv(1)      ! Maximum number of steps

      read(55,*) dt_variable    ! Variable or fixed dt
      if (dt_variable) then
         method(1) = 1
      else
         method(1) = 0
      end if
      read(55,*) method(2)    ! Order
      ! method(3) (transverse order) not used in 1D
      ! No dimensional splitting in 1D
      read(55,*) method(4)    ! Verbosity
      read(55,*) method(5)    ! Source term splitting style
      read(55,*) method(6)    ! Index into aux of capacity function
      method(7) = maux        ! Number of aux variables

      read(55,*) use_fwaves

      allocate(mthlim(mwaves), stat=allocate_status)
      if (allocate_status .ne. 0) then
         print *, '*** Error allocating mthlim array; exiting claw1ez'
         go to 900
      end if
      read(55,*) (mthlim(i), i = 1, mwaves)

      read(55,*) mbc
      read(55,*) mthbc(1)
      read(55,*) mthbc(2)

      ! No restart in 1D, and there's nothing after the restart info, so
      ! just close the file now.
      close(unit=55)


      ! Check consistency for periodic BCs
      if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) ' require mthbc(1) and mthbc(2) BOTH be set to 2'
         stop 
      endif


      ! Figure out size of work array needed
      mwork = (mx + 2*mbc) * (2 + 4*meqn + mwaves + meqn*mwaves)

c
c
      write(6,*) 'running...'
      write(6,*) ' '
c
c     # grid spacing
      dx = (xupper - xlower) / float(mx)
c

c     # time increments between outputing solution:
      if (outstyle .eq. 1) then
         dtout = (tfinal - t0)/float(nout)
         endif
c
c
c     # call user's routine setprob to set any specific parameters
c     # or other initialization required.
c
      call setprob

      ! Allocate aux
      if (maux > 0) then
         allocate(aux(maux, 1-mbc:mx+mbc), stat=allocate_status)
      else
         ! Allocate dummy array to prevent segfaults if it's
         ! dereferenced.  May be unnecessary.
         allocate(aux(1,1), stat=allocate_status)
      end if
      if (allocate_status .ne. 0) then
         print *, '*** Error allocating aux array; exiting claw1ez'
         go to 900
      end if

c        
c     # set aux array:
c
      if (maux .gt. 0)  then
         call setaux(mbc,mx,xlower,dx,maux,aux)
      endif

      ! Allocate q
      allocate(q(meqn, 1-mbc:mx+mbc), stat=allocate_status)
      if (allocate_status .ne. 0) then
         print *, '*** Error allocating q array; exiting claw1ez'
         go to 900
      end if

c
c     # set initial conditions:
c
      call qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c
c        # output initial data
      call out1(meqn,mbc,mx,xlower,dx,q,t0,0,aux,maux,
     &          outaux_init_only .or. outaux_always)

      ! Allocate work array
      allocate(work(mwork), stat=allocate_status)
      if (allocate_status .ne. 0) then
         print *, '*** Error allocating work array; exiting claw1ez'
         go to 900
      end if
c
c     ----------
c     Main loop:
c     ----------
c
      tend = t0
      do n=1,nout
         tstart = tend
         if (outstyle .eq. 1)  tend = tstart + dtout
         if (outstyle .eq. 2)  tend = tout(n)
         if (outstyle .eq. 3)  tend = tstart - 1.d0  !# single-step mode
c
         call claw1(meqn,mwaves,maux,mbc,mx,
     &           q,aux,xlower,dx,tstart,tend,dtv,cflv,nv,method,mthlim,
     &           mthbc,work,mwork,use_fwaves,info,bc1,rp1,src1,b4step1)
c
c        # check to see if an error occured:
         if (info .ne. 0) then
            write(6,*) '*** ERROR in claw1 ***  info =',info
            if (info.eq.1) then
               write(6,*) '***   either mx > maxmx or mbc < 2'
               endif
            if (info.eq.2) then
               write(6,*) '***   dt does not divide (tend - tstart)'
               write(6,*) '***   and dt is fixed since method(1)=0'
               endif
            if (info.eq.3) then
               write(6,*) '***   method(1)=1 and cflv(2) > cflv(1)'
               endif
            if (info.eq.4) then
               write(6,*) '***   mwork is too small'
               endif
            if (info.eq.11) then
               write(6,*) '***   Too many times steps, n > nv(1)'
               endif
            if (info.eq.12) then
               write(6,*) 
     &          '***   The Courant number is greater than cflv(1)'
               write(6,*) '***   and dt is fixed since method(1)=0'
               endif

            go to 900
            endif
c
         dtv(1) = dtv(5)  !# use final dt as starting value on next call
c
c        # output solution at this time
c        ------------------------------
c
c        # if outstyle=1 or 2, then nstepout=1 and we output every time
c        # we reach this point, since claw1 was called for the entire time
c        # increment between outputs.
c
c        # if outstyle=3 then we only output if we have taken nstepout
c        # time steps since the last output.

c        # iframe is the frame number used to form file names in out1
         iframe = n/nstepout
	 if (iframe*nstepout .eq. n) then
            call out1(meqn,mbc,mx,xlower,dx,q,tend,iframe,
     &                aux,maux,outaux_always)
            write(6,601) iframe,tend
            write(10,1010) tend,info,dtv(3),dtv(4),dtv(5),
     &           cflv(3),cflv(4),nv(2)
	    endif
c
c        # formats for writing out information about this call to claw:
c
  601    format('CLAW1EZ: Frame ',i4,
     &           ' output files done at time t =',
     &           d12.4,/)
c
 1010    format('tend =',d15.4,/,
     &       'info =',i5,/,'smallest dt =',d15.4,/,'largest dt =',
     &       d15.4,/,'last dt =',d15.4,/,'largest cfl =',
     &         d15.4,/,'last cfl =',d15.4,/,'steps taken =',i4,/)
c
      end do
c
  900 continue
      if (allocated(q))        deallocate(q)
      if (allocated(aux))      deallocate(aux)
      if (allocated(work))     deallocate(work)
      if (allocated(mthlim))   deallocate(mthlim)
      if (allocated(tout))     deallocate(tout)
      if (allocated(iout_q))   deallocate(iout_q)
      if (allocated(iout_aux)) deallocate(iout_aux)

c
      return 
      end

