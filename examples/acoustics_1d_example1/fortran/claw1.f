c
c     ==============================================================
      subroutine claw1(meqn,mwaves,maux,mbc,mx,
     &           q,aux,xlower,dx,tstart,tend,dtv,cflv,nv,method,mthlim,
     &           mthbc,work,mwork,use_fwave,info,bc1,rp1,src1,b4step1)
c     ==============================================================
c
c  Solves a hyperbolic system of conservation laws in one space dimension
c  of the general form 
c
c     capa * q_t + A q_x = psi
c
c  The "capacity function" capa(x) and source term psi are optional
c  (see below).
c
c  For a more complete description see the documentation at
c      http://www.amath.washington.edu/~claw
c
c  Sample driver programs and user-supplied subroutines are available.
c  See the the directories claw/clawpack/1d/example* for some examples, and
c  codes in claw/applications for more extensive examples.
c
c  --------------------------------------------------------
c
c  The user must supply the following subroutines:
c
c    bc1, rp1        subroutines specifying the boundary conditions and 
c                    Riemann solver.
c                    These are described in greater detail below.
c
c    b4step1            The routine b4step1 is called each time step and
c                       can be supplied by the user in order to perform
c                       other operations that are necessary every time
c                       step.  For example, if the variables stored in
c                       the aux arrays are time-dependent then these
c                       values can be set.
c
c  In addition, if the equation contains source terms psi, then the user
c  must provide:
c
c    src1               subroutine that solves capa * q_t = psi
c                       over a single time step.
c
c  These routines must be declared EXTERNAL in the main program.
c  For description of the calling sequences, see below.
c
c  Dummy routines b4step1.f and src1.f are available in 
c       claw/clawpack/1d/lib
c
c
c
c  Description of parameters...
c  ----------------------------
c
c
c    meqn is the number of equations in the system of
c         conservation laws.
c
c    mwaves is the number of waves that result from the
c           solution of each Riemann problem.  Often mwaves = meqn but
c           for some problems these may be different.
c
c    mbc is the number of "ghost cells" that must be added on to each
c       side of the domain to handle boundary conditions.  The cells
c       actually in the physical domain are labelled from 1 to mx in x.
c       The arrays are dimensioned actually indexed from 1-mbc to mx+mbc.
c       For the methods currently implemented, mbc = 2 should be used.
c       If the user implements another method that has a larger stencil and
c       hence requires more ghost cells, a larger value of mbc could be used.
c       q is extended from the physical domain to the ghost cells by the
c       user-supplied routine bc1.
c
c    mx is the number of grid cells in the x-direction, in the
c       physical domain.  In addition there are mbc grid cells
c       along each edge of the grid that are used for boundary
c       conditions.
c 
c    q(meqn, 1-mbc:mx+mbc) 
c        On input:  initial data at time tstart.
c        On output: final solution at time tend.
c        q(m,i) = value of mth component in the i'th cell.
c        Values within the physical domain are in q(m,i) 
c                for i = 1,2,...,mx
c        mbc extra cells on each end are needed for boundary conditions
c        as specified in the routine bc1.
c
c    aux(maux, 1-mbc:mx+mbc)
c        Array of auxiliary variables that are used in specifying the problem.
c        If method(7) = 0 then there are no auxiliary variables and aux
c                         can be a dummy variable.
c        If method(7) = maux > 0 then there are maux auxiliary variables
c                         and aux must be dimensioned as above.
c
c        Capacity functions are one particular form of auxiliary variable.
c        These arise in some applications, e.g. variable coefficients in
c        advection or acoustics problems.
c        See Clawpack Note # 5 for examples.
c
c        If method(6) = 0 then there is no capacity function.
c        If method(6) = mcapa > 0  then there is a capacity function and
c            capa(i), the "capacity" of the i'th cell, is assumed to be
c            stored in aux(mcapa,i).
c            In this case we require method(7).ge.mcapa.
c
c    dx = grid spacing in x.  
c         (for a computation in ax <= x <= bx,  set dx = (bx-ax)/mx.)
c
c    tstart = initial time.
c
c    tend = Desired final time (on input).
c              If tend<tstart, then claw1 returns after a single successful
c                 time step has been taken (single-step mode).
c              Otherwise, as many steps are taken as needed to reach tend, 
c                 up to a maximum of nv(1).
c         = Actual time reached (on output).
c
c    dtv(1:5) = array of values related to the time step:
c               (Note: method(1)=1 indicates variable size time steps)
c         dtv(1) = value of dt to be used in all steps if method(1) = 0
c                = value of dt to use in first step if method(1) = 1
c         dtv(2) = unused if method(1) = 0.
c                = maximum dt allowed if method(1) = 1.
c         dtv(3) = smallest dt used (on output)
c         dtv(4) = largest dt used (on output)
c         dtv(5) = dt used in last step (on output)
c
c    cflv(1:4) = array of values related to Courant number:
c         cflv(1) = maximum Courant number to be allowed.  With variable
c                   time steps the step is repeated if the Courant
c                   number is larger than this value.  With fixed time
c                   steps the routine aborts.  Usually cflv(1)=1.0
c                   should work.
c         cflv(2) = unused if method(1) = 0.
c                 = desired Courant number if method(1) = 1.
c                   Should be somewhat less than cflv(1), e.g. 0.9
c         cflv(3) = largest Courant number observed (on output).
c         cflv(4) = Courant number in last step (on output).
c
c    nv(1:2) = array of values related to the number of time steps:
c         nv(1) = unused if method(1) = 0
c               = maximum number of time steps allowed if method(1) = 1
c         nv(2) = number of time steps taken (on output).
c
c    method(1:7) = array of values specifying the numerical method to use
c         method(1) = 0 if fixed size time steps are to be taken.
c                       In this case, dt = dtv(1) in all steps.
c                   = 1 if variable time steps are to be used.
c                       In this case, dt = dtv(1) in the first step and
c                       thereafter the value cflv(2) is used to choose the
c                       next time step based on the maximum wave speed seen
c                       in the previous step.  Note that since this value
c                       comes from the previous step, the Courant number will
c                       not in general be exactly equal to the desired value
c                       If the actual Courant number in the next step is
c                       greater than 1, then this step is redone with a 
c                       smaller dt.
c
c         method(2) = 1 if Godunov's method is to be used, with no 2nd order
c                       corrections.
c                   = 2 if second order correction terms are to be added, with
c                       a flux limiter as specified by mthlim.  
c
c         method(3)  is not used in one-dimension.
c
c         method(4) = 0 to suppress printing
c                   = 1 to print dt and Courant number every time step
c
c         method(5) = 0 if there is no source term psi.  In this case
c                       the subroutine src1 is never called so a dummy
c                       parameter can be given.
c                   = 1 if there is a source term.  In this case
c                       the subroutine src1 must be provided and a
c                       fractional step method is used.
c                       In each time step the following sequence is followed:
c                            call bc to extend data to ghost cells
c                            call step1 to advance hyperbolic eqn by dt
c                            call src1 to advance source terms by dt
c                   = 2 if there is a source term and Strang splitting is to
c                       be used instead of the Godunov splitting above.
c                       In each time step the following sequence is followed:
c                            call bc to extend data to ghost cells
c                            call src1 to advance source terms by dt/2
c                            call step1 to advance hyperbolic equation by dt
c                            call src1 to advance source terms by dt/2
c                       For most problems 1 is recommended rather than 2
c                       since it is less expensive and works essentially as
c                       well on most problems.

c
c         method(6) = 0 if there is no capacity function capa.
c                   = mcapa > 0 if there is a capacity function.  In this case
c                       aux(i,mcapa) is the capacity of the i'th cell and you
c                       must also specify method(7) .ge. mcapa and set aux.
c
c         method(7) = 0 if there is no aux array used.
c                   = maux > 0  if there are maux auxiliary variables.
c
c
c         The recommended choice of methods for most problems is
c            method(1) = 1,  method(2) = 2.
c
c    mthlim(1:mwaves) = array of values specifying the flux limiter to be used
c                     in each wave family mw.  Often the same value will be used
c                     for each value of mw, but in some cases it may be
c                     desirable to use different limiters.  For example,
c                     for the Euler equations the superbee limiter might be
c                     used for the contact discontinuity (mw=2) while another
c                     limiter is used for the nonlinear waves.  Several limiters
c                     are built in and others can be added by modifying the
c                     subroutine philim.
c
c        mthlim(mw) = 0 for no limiter
c                   = 1 for minmod
c                   = 2 for superbee
c                   = 3 for van Leer
c                   = 4 for monotonized centered
c
c
c    work(mwork) = double precision work array of length at least mwork
c
c    mwork = length of work array.  Must be at least
c               (mx + 2*mbc) * (2 + 4*meqn + mwaves + meqn*mwaves)
c            If mwork is too small then the program returns with info = 4
c            and prints the necessary value of mwork to unit 6.
c
c            
c    info = output value yielding error information:
c         = 0 if normal return.
c         = 1 if mbc.lt.2
c         = 2 if method(1)=0 and dt doesn't divide (tend - tstart).
c         = 3 if method(1)=1 and cflv(2) > cflv(1).
c         = 4 if mwork is too small.
c         = 11 if the code attempted to take too many time steps, n > nv(1).
c              This could only happen if method(1) = 1 (variable time steps).
c         = 12 if the method(1)=0 and the Courant number is greater than 1
c              in some time step.
c
c           Note: if info.ne.0, then tend is reset to the value of t actually
c           reached and q contains the value of the solution at this time.
c
c    User-supplied subroutines
c    -------------------------
c
c    bc1 = subroutine that specifies the boundary conditions.  
c         This subroutine should extend the values of q from cells
c         1:mx to the mbc ghost cells along each edge of the domain.
c
c          The form of this subroutine is
c  -------------------------------------------------
c     subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
c     implicit double precision (a-h,o-z)
c     dimension   q(meqn, 1-mbc:mx+mbc)
c     dimension aux(maux, 1-mbc:mx+mbc)
c     dimension mthbc(2)
c  -------------------------------------------------
c
c    The routine claw/clawpack/1d/lib/bc1.f can be used to specify
c    various standard boundary conditions.
c
c
c    rp1 = user-supplied subroutine that implements the Riemann solver
c
c          The form of this subroutine is
c  -------------------------------------------------
c     subroutine rp1(meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
c     implicit double precision (a-h,o-z)
c     dimension   ql(meqn, 1-mbc:mx+mbc)
c     dimension   qr(meqn, 1-mbc:mx+mbc)
c     dimension auxl(maux, 1-mbc:mx+mbc)
c     dimension auxr(maux, 1-mbc:mx+mbc)
c     dimension wave(meqn, mwaves, 1-mbc:mx+mbc)
c     dimension    s(mwaves, 1-mbc:mx+mbc)
c     dimension amdq(meqn, 1-mbc:mx+mbc)
c     dimension apdq(meqn, 1-mbc:mx+mbc)
c  -------------------------------------------------
c
c         On input, ql contains the state vector at the left edge of each cell
c                   qr contains the state vector at the right edge of each cell
c                 auxl contains auxiliary values at the left edge of each cell
c                 auxr contains auxiliary values at the right edge of each cell
c
c         Note that the i'th Riemann problem has left state qr(:,i-1)
c                                            and right state ql(:,i)
c         In the standard clawpack routines, this Riemann solver is
c         called with ql=qr=q along this slice.  More flexibility is allowed
c         in case the user wishes to implement another solution method
c         that requires left and rate states at each interface.

c         If method(7)=maux > 0 then the auxiliary variables along this slice
c         are passed in using auxl and auxr.  Again, in the standard routines
c         auxl=auxr=aux in the call to rp1.
c
c          On output, 
c              wave(m,mw,i) is the m'th component of the jump across
c                              wave number mw in the ith Riemann problem.
c              s(mw,i) is the wave speed of wave number mw in the
c                              ith Riemann problem.
c              amdq(m,i) = m'th component of A^- Delta q,
c              apdq(m,i) = m'th component of A^+ Delta q,
c                     the decomposition of the flux difference
c                         f(qr(i-1)) - f(ql(i))
c                     into leftgoing and rightgoing parts respectively.
c
c           It is assumed that each wave consists of a jump discontinuity
c           propagating at a single speed, as results, for example, from a
c           Roe approximate Riemann solver.  An entropy fix can be included
c           into the specification of amdq and apdq.
c
c    src1 = subroutine for the source terms that solves the equation
c               capa * q_t = psi 
c           over time dt.
c
c           If method(5)=0 then the equation does not contain a source
c           term and this routine is never called.  A dummy argument can
c           be used with many compilers, or provide a dummy subroutine that
c           does nothing (such a subroutine can be found in
c           claw/clawpack/1d/lib/src1.f)
c
c          The form of this subroutine is
c  -------------------------------------------------
c     subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
c     implicit double precision (a-h,o-z)
c     dimension   q(meqn, 1-mbc:mx+mbc)
c     dimension aux(maux, 1-mbc:mx+mbc)
c  -------------------------------------------------
c      If method(7)=0  or the auxiliary variables are not needed in this solver,
c      then the latter dimension statement can be omitted, but aux should
c      still appear in the argument list.
c
c      On input, q(m,i) contains the data for solving the
c                source term equation.
c      On output, q(m,i) should have been replaced by the solution to
c                 the source term equation after a step of length dt.
c
c
c      b4step1 = subroutine that is called from claw1 before each call to
c                step1.  Use to set time-dependent aux arrays or perform
c                other tasks which must be done every time step.
c
c          The form of this subroutine is
c      
c  -------------------------------------------------
c      subroutine b4step1(mbc,mx,meqn,q,xlower,dx,time,dt,maux,aux)
c      implicit double precision (a-h,o-z)
c      dimension   q(meqn, 1-mbc:mx+mbc)
c      dimension aux(maux, 1-mbc:mx+mbc)
c  -------------------------------------------------
c
c
c
c
c =========================================================================
c
c  Copyright 1994 -- 2002 R. J. LeVeque
c
c  This software is made available for research and instructional use only. 
c  You may copy and use this software without charge for these non-commercial
c  purposes, provided that the copyright notice and associated text is
c  reproduced on all copies.  For all other uses (including distribution of
c  modified versions), please contact the author at the address given below. 
c  
c  *** This software is made available "as is" without any assurance that it
c  *** will work for your purposes.  The software may in fact have defects, so
c  *** use the software at your own risk.
c
c  --------------------------------------
c    CLAWPACK Version 4.1,  August, 2002
c    Webpage: http://www.amath.washington.edu/~claw
c  --------------------------------------
c    Author:  Randall J. LeVeque
c             Applied Mathematics
c             Box 352420
c             University of Washington, 
c             Seattle, WA 98195-2420
c             rjl@amath.washington.edu
c =========================================================================
c
c
c
c            
c    ======================================================================
c    Beginning of claw1 code
c    ======================================================================
c 
      implicit double precision (a-h,o-z)
      external bc1,rp1,src1,b4step1
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      dimension work(mwork)
      dimension mthlim(mwaves),method(7),dtv(5),cflv(4),nv(2)
      dimension mthbc(2)
      common /comxt/ dtcom,dxcom,tcom
c
c
      info = 0
      t = tstart
      maxn = nv(1)
      dt = dtv(1)   !# initial dt
      cflmax = 0.d0
      dtmin = dt
      dtmax = dt
      nv(2) = 0
c
c     # check for errors in data:
c
c
      if (method(1) .eq. 0) then
c        # fixed size time steps.  Compute the number of steps:
         if (tend .lt. tstart) then
c             # single step mode
              maxn = 1
           else
              maxn = nint((tend - tstart)/dt)    ! Round to nearest int
              if (dabs(maxn*dt - (tend-tstart)) .gt.
     &                          1d-5*(tend-tstart)) then
c                # dt doesn't divide time interval integer number of times
                 info = 2
                 go to 900
                 endif
           endif
         endif
c
      if (method(1).eq.1 .and. cflv(2).gt.cflv(1)) then
         info = 3
         go to 900
         endif
c
c     # partition work array into pieces for passing into step1:
      i0f = 1
      i0wave = i0f + (mx + 2*mbc) * meqn
      i0s = i0wave + (mx + 2*mbc) * meqn * mwaves
      i0dtdx = i0s + (mx + 2*mbc) * mwaves
      i0qwork = i0dtdx + (mx + 2*mbc) 
      i0amdq = i0qwork + (mx + 2*mbc) * meqn
      i0apdq = i0amdq + (mx + 2*mbc) * meqn
      i0dtdx = i0apdq + (mx + 2*mbc) * meqn
      i0end = i0dtdx + (mx + 2*mbc) - 1
c
      if (mwork .lt. i0end) then
         write(6,*) 'mwork must be increased to ',i0end
         info = 4
         go to 900
         endif
c
c     -----------
c     # main loop
c     -----------
c
      if (maxn.eq.0) go to 900
      do 100 n=1,maxn
         told = t   !# time at beginning of time step.

c        # adjust dt to hit tend exactly if we're near end of computation
c        #  (unless tend < tstart, which is a flag to take only a single step)
         if (told+dt.gt.tend .and. tstart.lt.tend) dt = tend - told

         if (method(1).eq.1) then
c           # save old q in case we need to retake step with smaller dt:
            call copyq1(meqn,mbc,mx,q,work(i0qwork))
            endif
c           
   40    continue
         dt2 = dt / 2.d0
         thalf = t + dt2  !# midpoint in time for Strang splitting
         t = told + dt    !# time at end of step

c        # store dt and t in the common block comxt in case they are needed
c        # in the Riemann solvers (for variable coefficients)
         tcom = told
         dtcom = dt
         dxcom = dx
c
c        ------------------------------------------------------------------
c        # main steps in algorithm:
c        ------------------------------------------------------------------
c
c        # extend data from grid to bordering boundary cells:
         call bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,told,dt,mthbc)
c
c
c        # call user-supplied routine which might set aux arrays
c        # for this time step, for example.

         call b4step1(mbc,mx,meqn,q,
     &                xlower,dx,told,dt,maux,aux)
c
c
         if (method(5).eq.2) then
c            # with Strang splitting for source term:
             call src1(meqn,mbc,mx,xlower,dx,q,maux,aux,told,dt2)
             endif
c
c        # take a step on the homogeneous conservation law:
         call step1(meqn,mwaves,mbc,maux,mx,q,aux,dx,dt,
     &             method,mthlim,cfl,work(i0f),work(i0wave),
     &             work(i0s),work(i0amdq),work(i0apdq),work(i0dtdx),
     &             use_fwave,rp1)
c
         if (method(5).eq.2) then
c            # source terms over a second half time step for Strang splitting:
c            # Note it is not so clear what time t should be used here if
c            # the source terms are time-dependent!
             call src1(meqn,mbc,mx,xlower,dx,q,maux,aux,thalf,dt2)
             endif

         if (method(5).eq.1) then
c            # source terms over a full time step:
             call src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
             endif
c

c
c        ------------------------------------------------------------------
c
         if (method(4) .eq. 1) write(6,601) n,cfl,dt,t
  601    format('CLAW1... Step',i4,
     &                   '   Courant number =',f6.3,'  dt =',d12.4,
     &                   '  t =',d12.4)
c
         if (method(1) .eq. 1) then
c           # choose new time step if variable time step
            if (cfl .gt. 0.d0) then
                dt = dmin1(dtv(2), dt * cflv(2)/cfl)
                dtmin = dmin1(dt,dtmin)
                dtmax = dmax1(dt,dtmax)
              else
                dt = dtv(2)
              endif
            endif
c
c        # check to see if the Courant number was too large:
c
         if (cfl .le. cflv(1)) then
c               # accept this step
                cflmax = dmax1(cfl,cflmax)
              else
c               # reject this step
                t = told
                call copyq1(meqn,mbc,mx,work(i0qwork),q)
c
                if (method(4) .eq. 1) then
                   write(6,602) 
  602              format('CLAW1 rejecting step... ',
     &                         'Courant number too large')
                   endif
                if (method(1).eq.1) then
c                   # if variable dt, go back and take a smaller step
                    go to 40
                  else
c                   # if fixed dt, give up and return
                    cflmax = dmax1(cfl,cflmax)
                    go to 900
                  endif
               endif
c
c        # see if we are done:
         nv(2) = nv(2) + 1
         if (t .ge. tend) go to 900
c
  100    continue
c
  900  continue
c 
c      # return information
c
       if (method(1).eq.1 .and. t.lt.tend .and. nv(2) .eq. maxn) then
c         # too many timesteps
          info = 11
          endif
c
       if (method(1).eq.0 .and. cflmax .gt. cflv(1)) then
c         # Courant number too large with fixed dt
          info = 12
          endif
       tend = t
       cflv(3) = cflmax
       cflv(4) = cfl
       dtv(3) = dtmin
       dtv(4) = dtmax
       dtv(5) = dt
       return 
       end

