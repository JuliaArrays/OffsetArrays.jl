using OffsetArrays
#include("copyq1.jl")
include("step1.jl")

#
#     ================================================================================================================
function claw1(parms::CParam,meqn::Int, mwaves::Int, maux::Int, mbc::Int, mx::Int,
               q::OffsetArray{Float64,2}, aux::OffsetArray{Float64,2},
               xlower::Float64, dx::Float64, tstart::Float64, tend::Float64,
               dtv::Array{Float64,1}, cflv::Array{Float64,1},
               nv::Array{Int,1}, method::Array{Int,1}, mthlim::Array{Int,1}, mthbc::Array{Int,1},
               #work,mwork, # ???
               use_fwave::Bool,
               # higher-order functions
               bc1, rp1, src1, b4step1)

#     ================================================================================================================
#
#  Solves a hyperbolic system of conservation laws in one space dimension
#  of the general form 
#
#     capa * q_t + A q_x = psi
#
#  The "capacity function" capa(x) and source term psi are optional
#  (see below).
#
#  For a more complete description see the documentation at
#      http://www.amath.washington.edu/~claw
#
#  Sample driver programs and user-supplied subroutines are available.
#  See the the directories claw/clawpack/1d/example* for some examples, and
#  codes in claw/applications for more extensive examples.
#
#  --------------------------------------------------------
#
#  The user must supply the following subroutines:
#
#    bc1, rp1        subroutines specifying the boundary conditions and 
#                    Riemann solver.
#                    These are described in greater detail below.
#
#    b4step1            The routine b4step1 is called each time step and
#                       can be supplied by the user in order to perform
#                       other operations that are necessary every time
#                       step.  For example, if the variables stored in
#                       the aux arrays are time-dependent then these
#                       values can be set.
#
#  In addition, if the equation contains source terms psi, then the user
#  must provide:
#
#    src1               subroutine that solves capa * q_t = psi
#                       over a single time step.
#
#  These routines must be declared EXTERNAL in the main program.
#  For description of the calling sequences, see below.
#
#  Dummy routines b4step1.f and src1.f are available in 
#       claw/clawpack/1d/lib
#
#
#
#  Description of parameters...
#  ----------------------------
#
#
#    meqn is the number of equations in the system of
#         conservation laws.
#
#    mwaves is the number of waves that result from the
#           solution of each Riemann problem.  Often mwaves = meqn but
#           for some problems these may be different.
#
#    mbc is the number of "ghost cells" that must be added on to each
#       side of the domain to handle boundary conditions.  The cells
#       actually in the physical domain are labelled from 1 to mx in x.
#       The arrays are dimensioned actually indexed from 1-mbc to mx+mbc.
#       For the methods currently implemented, mbc = 2 should be used.
#       If the user implements another method that has a larger stencil and
#       hence requires more ghost cells, a larger value of mbc could be used.
#       q is extended from the physical domain to the ghost cells by the
#       user-supplied routine bc1.
#
#    mx is the number of grid cells in the x-direction, in the
#       physical domain.  In addition there are mbc grid cells
#       along each edge of the grid that are used for boundary
#       conditions.
# 
#    q(meqn, 1-mbc:mx+mbc) 
#        On input:  initial data at time tstart.
#        On output: final solution at time tend.
#        q(m,i) = value of mth component in the i'th cell.
#        Values within the physical domain are in q(m,i) 
#                for i = 1,2,...,mx
#        mbc extra cells on each end are needed for boundary conditions
#        as specified in the routine bc1.
#
#    aux(maux, 1-mbc:mx+mbc)
#        Array of auxiliary variables that are used in specifying the problem.
#        If method(7) = 0 then there are no auxiliary variables and aux
#                         can be a dummy variable.
#        If method(7) = maux > 0 then there are maux auxiliary variables
#                         and aux must be dimensioned as above.
#
#        Capacity functions are one particular form of auxiliary variable.
#        These arise in some applications, e.g. variable coefficients in
#        advection or acoustics problems.
#        See Clawpack Note # 5 for examples.
#
#        If method(6) = 0 then there is no capacity function.
#        If method(6) = mcapa > 0  then there is a capacity function and
#            capa(i), the "capacity" of the i'th cell, is assumed to be
#            stored in aux(mcapa,i).
#            In this case we require method(7).ge.mcapa.
#
#    dx = grid spacing in x.  
#         (for a computation in ax <= x <= bx,  set dx = (bx-ax)/mx.)
#
#    tstart = initial time.
#
#    tend = Desired final time (on input).
#              If tend<tstart, then claw1 returns after a single successful
#                 time step has been taken (single-step mode).
#              Otherwise, as many steps are taken as needed to reach tend, 
#                 up to a maximum of nv(1).
#         = Actual time reached (on output).
#
#    dtv(1:5) = array of values related to the time step:
#               (Note: method(1)=1 indicates variable size time steps)
#         dtv(1) = value of dt to be used in all steps if method(1) = 0
#                = value of dt to use in first step if method(1) = 1
#         dtv(2) = unused if method(1) = 0.
#                = maximum dt allowed if method(1) = 1.
#         dtv(3) = smallest dt used (on output)
#         dtv(4) = largest dt used (on output)
#         dtv(5) = dt used in last step (on output)
#
#    cflv(1:4) = array of values related to Courant number:
#         cflv(1) = maximum Courant number to be allowed.  With variable
#                   time steps the step is repeated if the Courant
#                   number is larger than this value.  With fixed time
#                   steps the routine aborts.  Usually cflv(1)=1.0
#                   should work.
#         cflv(2) = unused if method(1) = 0.
#                 = desired Courant number if method(1) = 1.
#                   Should be somewhat less than cflv(1), e.g. 0.9
#         cflv(3) = largest Courant number observed (on output).
#         cflv(4) = Courant number in last step (on output).
#
#    nv(1:2) = array of values related to the number of time steps:
#         nv(1) = unused if method(1) = 0
#               = maximum number of time steps allowed if method(1) = 1
#         nv(2) = number of time steps taken (on output).
#
#    method(1:7) = array of values specifying the numerical method to use
#         method(1) = 0 if fixed size time steps are to be taken.
#                       In this case, dt = dtv(1) in all steps.
#                   = 1 if variable time steps are to be used.
#                       In this case, dt = dtv(1) in the first step and
#                       thereafter the value cflv(2) is used to choose the
#                       next time step based on the maximum wave speed seen
#                       in the previous step.  Note that since this value
#                       comes from the previous step, the Courant number will
#                       not in general be exactly equal to the desired value
#                       If the actual Courant number in the next step is
#                       greater than 1, then this step is redone with a 
#                       smaller dt.
#
#         method(2) = 1 if Godunov's method is to be used, with no 2nd order
#                       corrections.
#                   = 2 if second order correction terms are to be added, with
#                       a flux limiter as specified by mthlim.  
#
#         method(3)  is not used in one-dimension.
#
#         method(4) = 0 to suppress printing
#                   = 1 to print dt and Courant number every time step
#
#         method(5) = 0 if there is no source term psi.  In this case
#                       the subroutine src1 is never called so a dummy
#                       parameter can be given.
#                   = 1 if there is a source term.  In this case
#                       the subroutine src1 must be provided and a
#                       fractional step method is used.
#                       In each time step the following sequence is followed:
#                            call bc to extend data to ghost cells
#                            call step1 to advance hyperbolic eqn by dt
#                            call src1 to advance source terms by dt
#                   = 2 if there is a source term and Strang splitting is to
#                       be used instead of the Godunov splitting above.
#                       In each time step the following sequence is followed:
#                            call bc to extend data to ghost cells
#                            call src1 to advance source terms by dt/2
#                            call step1 to advance hyperbolic equation by dt
#                            call src1 to advance source terms by dt/2
#                       For most problems 1 is recommended rather than 2
#                       since it is less expensive and works essentially as
#                       well on most problems.

#
#         method(6) = 0 if there is no capacity function capa.
#                   = mcapa > 0 if there is a capacity function.  In this case
#                       aux(i,mcapa) is the capacity of the i'th cell and you
#                       must also specify method(7) .ge. mcapa and set aux.
#
#         method(7) = 0 if there is no aux array used.
#                   = maux > 0  if there are maux auxiliary variables.
#
#
#         The recommended choice of methods for most problems is
#            method(1) = 1,  method(2) = 2.
#
#    mthlim(1:mwaves) = array of values specifying the flux limiter to be used
#                     in each wave family mw.  Often the same value will be used
#                     for each value of mw, but in some cases it may be
#                     desirable to use different limiters.  For example,
#                     for the Euler equations the superbee limiter might be
#                     used for the contact discontinuity (mw=2) while another
#                     limiter is used for the nonlinear waves.  Several limiters
#                     are built in and others can be added by modifying the
#                     subroutine philim.
#
#        mthlim(mw) = 0 for no limiter
#                   = 1 for minmod
#                   = 2 for superbee
#                   = 3 for van Leer
#                   = 4 for monotonized centered
#
#
#    work(mwork) = double precision work array of length at least mwork
#
#    mwork = length of work array.  Must be at least
#               (mx + 2*mbc) * (2 + 4*meqn + mwaves + meqn*mwaves)
#            If mwork is too small then the program returns with info = 4
#            and prints the necessary value of mwork to unit 6.
#
#            
#    info = output value yielding error information:
#         = 0 if normal return.
#         = 1 if mbc.lt.2
#         = 2 if method(1)=0 and dt doesn't divide (tend - tstart).
#         = 3 if method(1)=1 and cflv(2) > cflv(1).
#         = 4 if mwork is too small.
#         = 11 if the code attempted to take too many time steps, n > nv(1).
#              This could only happen if method(1) = 1 (variable time steps).
#         = 12 if the method(1)=0 and the Courant number is greater than 1
#              in some time step.
#
#           Note: if info.ne.0, then tend is reset to the value of t actually
#           reached and q contains the value of the solution at this time.
#
#    User-supplied subroutines
#    -------------------------
#
#    bc1 = subroutine that specifies the boundary conditions.  
#         This subroutine should extend the values of q from cells
#         1:mx to the mbc ghost cells along each edge of the domain.
#
#          The form of this subroutine is
#  -------------------------------------------------
#     subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
#     implicit double precision (a-h,o-z)
#     dimension   q(meqn, 1-mbc:mx+mbc)
#     dimension aux(maux, 1-mbc:mx+mbc)
#     dimension mthbc(2)
#  -------------------------------------------------
#
#    The routine claw/clawpack/1d/lib/bc1.f can be used to specify
#    various standard boundary conditions.
#
#
#    rp1 = user-supplied subroutine that implements the Riemann solver
#
#          The form of this subroutine is
#  -------------------------------------------------
#     subroutine rp1(meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
#     implicit double precision (a-h,o-z)
#     dimension   ql(meqn, 1-mbc:mx+mbc)
#     dimension   qr(meqn, 1-mbc:mx+mbc)
#     dimension auxl(maux, 1-mbc:mx+mbc)
#     dimension auxr(maux, 1-mbc:mx+mbc)
#     dimension wave(meqn, mwaves, 1-mbc:mx+mbc)
#     dimension    s(mwaves, 1-mbc:mx+mbc)
#     dimension amdq(meqn, 1-mbc:mx+mbc)
#     dimension apdq(meqn, 1-mbc:mx+mbc)
#  -------------------------------------------------
#
#         On input, ql contains the state vector at the left edge of each cell
#                   qr contains the state vector at the right edge of each cell
#                 auxl contains auxiliary values at the left edge of each cell
#                 auxr contains auxiliary values at the right edge of each cell
#
#         Note that the i'th Riemann problem has left state qr(:,i-1)
#                                            and right state ql(:,i)
#         In the standard clawpack routines, this Riemann solver is
#         called with ql=qr=q along this slice.  More flexibility is allowed
#         in case the user wishes to implement another solution method
#         that requires left and rate states at each interface.

#         If method(7)=maux > 0 then the auxiliary variables along this slice
#         are passed in using auxl and auxr.  Again, in the standard routines
#         auxl=auxr=aux in the call to rp1.
#
#          On output, 
#              wave(m,mw,i) is the m'th component of the jump across
#                              wave number mw in the ith Riemann problem.
#              s(mw,i) is the wave speed of wave number mw in the
#                              ith Riemann problem.
#              amdq(m,i) = m'th component of A^- Delta q,
#              apdq(m,i) = m'th component of A^+ Delta q,
#                     the decomposition of the flux difference
#                         f(qr(i-1)) - f(ql(i))
#                     into leftgoing and rightgoing parts respectively.
#
#           It is assumed that each wave consists of a jump discontinuity
#           propagating at a single speed, as results, for example, from a
#           Roe approximate Riemann solver.  An entropy fix can be included
#           into the specification of amdq and apdq.
#
#    src1 = subroutine for the source terms that solves the equation
#               capa * q_t = psi 
#           over time dt.
#
#           If method(5)=0 then the equation does not contain a source
#           term and this routine is never called.  A dummy argument can
#           be used with many compilers, or provide a dummy subroutine that
#           does nothing (such a subroutine can be found in
#           claw/clawpack/1d/lib/src1.f)
#
#          The form of this subroutine is
#  -------------------------------------------------
#     subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
#     implicit double precision (a-h,o-z)
#     dimension   q(meqn, 1-mbc:mx+mbc)
#     dimension aux(maux, 1-mbc:mx+mbc)
#  -------------------------------------------------
#      If method(7)=0  or the auxiliary variables are not needed in this solver,
#      then the latter dimension statement can be omitted, but aux should
#      still appear in the argument list.
#
#      On input, q(m,i) contains the data for solving the
#                source term equation.
#      On output, q(m,i) should have been replaced by the solution to
#                 the source term equation after a step of length dt.
#
#
#      b4step1 = subroutine that is called from claw1 before each call to
#                step1.  Use to set time-dependent aux arrays or perform
#                other tasks which must be done every time step.
#
#          The form of this subroutine is
#      
#  -------------------------------------------------
#      subroutine b4step1(mbc,mx,meqn,q,xlower,dx,time,dt,maux,aux)
#      implicit double precision (a-h,o-z)
#      dimension   q(meqn, 1-mbc:mx+mbc)
#      dimension aux(maux, 1-mbc:mx+mbc)
#  -------------------------------------------------
#
#
#
#
# =========================================================================
#
#  Copyright 1994 -- 2002 R. J. LeVeque
#
#  This software is made available for research and instructional use only. 
#  You may copy and use this software without charge for these non-commercial
#  purposes, provided that the copyright notice and associated text is
#  reproduced on all copies.  For all other uses (including distribution of
#  modified versions), please contact the author at the address given below. 
#  
#  *** This software is made available "as is" without any assurance that it
#  *** will work for your purposes.  The software may in fact have defects, so
#  *** use the software at your own risk.
#
#  --------------------------------------
#    CLAWPACK Version 4.1,  August, 2002
#    Webpage: http://www.amath.washington.edu/~claw
#  --------------------------------------
#    Author:  Randall J. LeVeque
#             Applied Mathematics
#             Box 352420
#             University of Washington, 
#             Seattle, WA 98195-2420
#             rjl@amath.washington.edu
# =========================================================================
#
#
#
#            
#    ======================================================================
#    Beginning of claw1 code
#    ======================================================================
# 
#      implicit double precision (a-h,o-z)
#      external bc1,rp1,src1,b4step1
#      dimension q(meqn,1-mbc:mx+mbc)
#      dimension aux(maux,1-mbc:mx+mbc)
#      dimension work(mwork)
#      dimension mthlim(mwaves),method(7),dtv(5),cflv(4),nv(2)
#      dimension mthbc(2)
#      common /comxt/ dtcom,dxcom,tcom
#
#
    
    info::Int = 0
    q_save = OffsetArray(Float64, 1:meqn, 1-mbc:mx+mbc)

    info = 0
    t = tstart
    maxn = nv[1]
    dt = dtv[1]   # initial dt
    cflmax = 0.0
    dtmin = dt
    dtmax = dt
    nv[2] = 0

    # check for errors in data:

    if method[1] == 0
        # fixed size time steps.  Compute the number of steps:
        if tend < tstart
            # single step mode
            maxn = 1
        else
            maxn = convert(Uint, (tend - tstart)/dt)    # Round to nearest int
            if abs(maxn*dt - (tend-tstart)) >
                        1e-5*(tend-tstart)
                # dt doesn't divide time interval integer number of times
                info = 2
                # [local goto #101](https://github.com/JuliaLang/julia/issues/101)
                @goto l900
            end
        end
    end

    if method[1] == 1 && cflv[2] > cflv[1]
        info = 3
        @goto l900
    end

    #
    #     -----------
    #     # main loop
    #     -----------
    #

    if maxn == 0
        return
    end

    f    = OffsetArray(Float64, 1:meqn,              1-mbc:mx+mbc)
    s    = OffsetArray(Float64,            1:mwaves, 1-mbc:mx+mbc)
    wave = OffsetArray(Float64, 1:meqn,    1:mwaves, 1-mbc:mx+mbc)
    amdq = OffsetArray(Float64, 1:meqn,              1-mbc:mx+mbc)
    apdq = OffsetArray(Float64, 1:meqn,              1-mbc:mx+mbc)
    dtdx = OffsetArray(Float64,                      1-mbc:mx+mbc)

    cfl::Float64 = 0.0
    for n=1:maxn
        told = t   # time at beginning of time step.

        # adjust dt to hit tend exactly if we're near end of computation
        #  (unless tend < tstart, which is a flag to take only a single step)
        if told+dt > tend && tstart < tend
            dt = tend - told
        end

        if method[1] == 1
            # save old q in case we need to retake step with smaller dt:
            #copyq1(meqn,mbc,mx,q,work(i0qwork))
            q_save = q
        end

    #           
    @label l40

        dt2 = dt / 2.0
        #println("dt2: $dt2")
        thalf = t + dt2  # midpoint in time for Strang splitting
        t = told + dt    # time at end of step

        # store dt and t in the common block comxt in case they are needed
        # in the Riemann solvers (for variable coefficients)
        tcom = told
        dtcom = dt
        dxcom = dx

        # ------------------------------------------------------------------
        # main steps in algorithm:
        # ------------------------------------------------------------------

        # extend data from grid to bordering boundary cells:
        bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,told,dt,mthbc)

        # call user-supplied routine which might set aux arrays
        # for this time step, for example.

        b4step1(mbc,mx,meqn,q,xlower,dx,told,dt,maux,aux)

        if method[5] == 2
            # with Strang splitting for source term:
            src1(meqn,mbc,mx,xlower,dx,q,maux,aux,told,dt2)
        end

        # take a step on the homogeneous conservation law:

        cfl =
        step1(parms,meqn,mwaves,mbc,maux,mx,q,aux,dx,dt,method,mthlim,
              f,wave,s,amdq,apdq,dtdx,
              use_fwave,rp1)

        if method[5] == 2
            # source terms over a second half time step for Strang splitting:
            # Note it is not so clear what time t should be used here if
            # the source terms are time-dependent!
            src1(meqn,mbc,mx,xlower,dx,q,maux,aux,thalf,dt2)
        end

        if method[5] == 1
            # source terms over a full time step:
            src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
        end

        #
        #        ------------------------------------------------------------------
        #
        if method[4] == 1
            @printf("CLAW1... Step %4d\n   Courant number =',%6.3f\n  dt =%12.4e,\n  t =',%12.4e\n",n,cfl,dt,t)
        end

        if method[1] == 1
            # choose new time step if variable time step
            if cfl > 0.0
                dt = min(dtv[2], dt * cflv[2]/cfl)
                dtmin = min(dt,dtmin)
                dtmax = max(dt,dtmax)
            else
                dt = dtv[2]
            end
        end

        # check to see if the Courant number was too large:

        if cfl <= cflv[1]
            # accept this step
            cflmax = max(cfl,cflmax)
            # as the step accepted exit the loop
            #println("ACCEPT: cfl : $cfl cflv[1] : $(cflv[1])")
        else
            #println("REJECT: cfl : $cfl cflv[1] : $(cflv[1])")
            # reject this step
            t = told
            #copyq1(meqn,mbc,mx,work(i0qwork),q)
            q = q_save

            if method[4] == 1
                println("CLAW1 rejecting step... ")
                println("Courant number too large")
            end

            if method[1] == 1
                # if variable dt, go back and take a smaller step
                @goto l40
            else # if fixed dt, give up and return
                cflmax = max(cfl,cflmax)
                break
            end

        end

        # see if we are done:
        nv[2] += 1

        println("t : $t tend : $tend")
        if t >= tend
            break
        end

    end


    #
@label l900
    
 
    # return information

    if method[1] == 1 && t < tend && nv[2] == maxn
        # too many timesteps
        info = 11
    end

    if method[1] == 0 && cflmax > cflv[1]
        # Courant number too large with fixed dt
        info = 12
    end
       
    tend = t
    cflv[3] = cflmax
    cflv[4] = cfl
    dtv[3] = dtmin
    dtv[4] = dtmax
    dtv[5] = dt

    return info 

end # claw1

