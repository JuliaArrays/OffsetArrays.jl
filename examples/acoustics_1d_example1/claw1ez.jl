
using OffsetArrays

include("ReadDataFile.jl")
include("SetProb.jl")

using ReadDataFile
using SetProb

include("qinit.jl")
include("out1.jl")
include("src1.jl")
include("bc1.jl")
include("rp1_acoustics.jl")
include("b4step1.jl")
include("claw1.jl")


#
#
#
#     ===============================================================
function claw1ez()    # No arguments
#     ===============================================================
#
#     An easy-to-use clawpack driver routine for simple applications
#     Documentation is available at
#                 http://www.amath.washington.edu/~claw/doc.html
#
#     Authors: Randall J. LeVeque, Grady I. Lemoine
#
#      implicit double precision (a-h,o-z)
#      external bc1,rp1,src1,b4step1

#      double precision, dimension(:,:), allocatable :: q, aux
#      double precision, dimension(:), allocatable :: work, tout
#      integer, dimension(:), allocatable :: mthlim, iout_q, iout_aux


#      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(2)
    method = Array(Int, 7)
    dtv    = Array(Float64, 5)
    cflv   = Array(Float64, 4)
    nv     = Array(Int, 2)
    mthbc  = Array(Int, 2)

###       integer :: allocate_status, outstyle
###       logical :: outaux_init_only, use_fwaves, output_t0
###       logical :: outaux_always, dt_variable

    outaux_init_only = false
    outaux_always = false

###       character*12 fname
### #
###       open(10,file='fort.info',status='unknown',form='formatted')
### 
### #
### #     # Read the input in standard form from claw.data:
### #     # For a description of input parameters see the documentation at
### #                 http://www.amath.washington.edu/~claw
### 
### #     ## New in 4.4:   Input file name changed from claw1ez.data
### #     ## Open file and skip over leading lines with # comments:
###       fname = 'claw.data'
###       call opendatafile(55,fname)
### 
### #
### #     # Read the input in standard form from claw.data:
### #
### 
###       read(55,*) ndim
### 
###       read(55,*) xlower
###       read(55,*) xupper
###       read(55,*) mx
###       read(55,*) meqn
###       read(55,*) mwaves
###       read(55,*) maux
###       read(55,*) t0
### 
    r = Reader("claw.data")

    ndim     = readint(r)
    xlower   = readflt(r)
    xupper   = readflt(r)
    mx       = readint(r)
    meqn     = readint(r)
    mwaves   = readint(r)
    maux     = readint(r)
    t0       = readflt(r)

    println("read from claw.data : ndim: $ndim xlower: $xlower xupper: $xupper mx: $mx meqn: $meqn mwaves: $mwaves maux: $maux t0: $t0")

###       read(55,*) outstyle

    outstyle = readint(r)
 
###       if (outstyle == 1) then
###          read(55,*) nout
###          read(55,*) tfinal
###          read(55,*) output_t0    ! Not currently used
###          nstepout = 1

    tout = Array(Float64, 1)
    if     outstyle == 1
        nout      = readint(r)
        tfinal    = readflt(r)
        output_t0 = readbool(r)       # Not currently used
        nstepout  = 1

        println("nout: $nout tfinal: $tfinal output_t0: $output_t0")


    elseif outstyle == 2
        nout      = readint(r)
        tout      = readfltarray(r, nout)
        println("tout: $tout")
        nstepout  = 1


###       else if (outstyle == 3) then
###          read(55,*) nstepout
###          read(55,*) nstop
###          read(55,*) output_t0
###          nout = nstop
###       else
###          print *, '*** Unrecognized output style ', outstyle
###          print *, '*** Exiting claw1ez'
###          go to 900
###       end if
### 

    elseif outstyle == 3

    else
        println("*** Unrecognized output style $outstyle")
        error("*** Exiting claw1ez")
    end

###       read(55,*) output_format    ! Not used yet
    output_format = readint(r)    # Not used yet

###       ! These iout variables are not currently used, but hang onto them
###       ! anyway in case somebody wants to use them at a future date.  The
###       ! same goes for outaux_init_only.
###       allocate(iout_q(meqn), stat=allocate_status)
###       if (allocate_status .ne. 0) then
###          print *, '*** Error allocating iout_q array; exiting claw1ez'
###          go to 900    ! Exception handling, old school style
###       end if
###       read(55,*) (iout_q(i), i = 1, meqn)
    inout_q = readfltarray(r, meqn)
    println("inout_q: $inout_q")

#    if maux > 0
###          allocate(iout_aux(maux), stat=allocate_status)
###          if (allocate_status .ne. 0) then
###             print *, '*** Error allocating iout_aux array;',
###      &               ' exiting claw1ez'
###             go to 900
###          end if

#        iout_aux

###          read(55,*) (iout_aux(i), i = 1, maux)
###          read(55,*) outaux_init_only
###          ! Not implementing selective output of aux fields yet
###          if (any(iout_aux .ne. 0)) then
###             outaux_always = .not. outaux_init_only
###          else
###             outaux_always = .false.
###             outaux_init_only = .false.
###          end if
###       else
###          outaux_always = .false.
###          outaux_init_only = .false.    ! Just to initialize
#    end


    iout_aux = Array(Int, 1)
    if maux > 0
        iout_aux = readintarray(r, maux)
        outaux_init_only = readbool(r)

        # TODO
    else
        outaux_always = false
        outaux_init_only = false    # Just to initialize
    end

    dtv[1]  = readflt(r)    # Initial dt
    dtv[2]  = readflt(r)    # Max dt
    cflv[1] = readflt(r)    # Max CFL number
    cflv[2] = readflt(r)    # Desired CFL number
    nv[1]   = readint(r)    # Maximum number of steps

    println("dtv: $dtv cflv: $cflv nv: $nv")

###       read(55,*) dt_variable    ! Variable or fixed dt
###       if (dt_variable) then
###          method(1) = 1
###       else
###          method(1) = 0
###       end if

    dt_variable::Bool = readbool(r)
    method[1] = dt_variable ? 1 : 0

    method[2] = readint(r) # Order

###       ! method(3) (transverse order) not used in 1D
###       ! No dimensional splitting in 1D
###       read(55,*) method(4)    ! Verbosity
###       read(55,*) method(5)    ! Source term splitting style
###       read(55,*) method(6)    ! Index into aux of capacity function
###       method(7) = maux        ! Number of aux variables
### 


    method[4] = readint(r) # Verbosity
    method[5] = readint(r) # Source term splitting style
    method[6] = readint(r) # Index into aux of capacity function
    method[7] = maux       # Number of aux variables

    println("method: $method")

###       read(55,*) use_fwaves
### 
    use_fwaves::Bool = readbool(r)

###       allocate(mthlim(mwaves), stat=allocate_status)
###       if (allocate_status .ne. 0) then
###          print *, '*** Error allocating mthlim array; exiting claw1ez'
###          go to 900
###       end if
###       read(55,*) (mthlim(i), i = 1, mwaves)

    mthlim::Array{Int, 1} = readintarray(r, mwaves)
    println("mthlim: $mthlim")

    mbc      = readint(r)
    mthbc[1] = readint(r)
    mthbc[2] = readint(r)

    println("mbc: $mbc mthbc: $mthbc")


    # No restart in 1D, and there's nothing after the restart info, so
    # just close the file now.
    close_reader(r) 

###       ! Check consistency for periodic BCs
###       if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
###      &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
###          write(6,*) '*** ERROR ***  periodic boundary conditions'
###          write(6,*) ' require mthbc(1) and mthbc(2) BOTH be set to 2'
###          stop 
###       endif

    # Check consistency for periodic BCs
    if (mthbc[1] == 2 && mthbc[2] != 2) ||
       (mthbc[2] == 2 && mthbc[1] != 2)
          error("*** ERROR ***  periodic boundary conditions require mthbc(1) and mthbc(2) BOTH be set to 2")
    end

###       ! Figure out size of work array needed
###       mwork = (mx + 2*mbc) * (2 + 4*meqn + mwaves + meqn*mwaves)
### 
### #
### #
###       write(6,*) 'running...'
###       write(6,*) ' '
### #
    println("running...")
    println(" ")

    # grid spacing
    dx = (xupper - xlower) / float(mx)

    # time increments between outputing solution:
    if outstyle == 1
        dtout = (tfinal - t0)/float(nout)
    end

    # call user's routine setprob to set any specific parameters
    # or other initialization required.
    #
    prob_data = CParam()
    read_prob_data(prob_data)

 
    # Allocate aux

    aux = (maux > 0) ? OffsetArray(Float64, 1:maux, 1-mbc:mx+mbc) : OffsetArray(Float64, 1:1,1:1)

### #        
### #     # set aux array:
### #
###       if (maux .gt. 0)  then
###          call setaux(mbc,mx,xlower,dx,maux,aux)
###       endif
### 
###       ! Allocate q
###       allocate(q(meqn, 1-mbc:mx+mbc), stat=allocate_status)
###       if (allocate_status .ne. 0) then
###          print *, '*** Error allocating q array; exiting claw1ez'
###          go to 900
###       end if

    q = OffsetArray(Float64, 1:meqn, 1-mbc:mx+mbc)

### #     # set initial conditions:
### #
###       call qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
### #

    #println("q before qinit call: $q")

    qinit(prob_data, meqn, mbc, mx, xlower, dx, q, maux, aux)

    #println("q after qinit call:")
    #Base.print(q)

    # output initial data
    out1(meqn,mbc,mx,xlower,dx,q,t0,0,aux,maux,
         (outaux_init_only || outaux_always) )

### 
###       ! Allocate work array
###       allocate(work(mwork), stat=allocate_status)
###       if (allocate_status .ne. 0) then
###          print *, '*** Error allocating work array; exiting claw1ez'
###          go to 900
###       end if


#
#     ----------
#     Main loop:
#     ----------
#
    tend = t0
    for n=1:nout
        tstart = tend
        if outstyle == 1
            tend = tstart + dtout
        elseif outstyle == 2
            tend = tout[n]
        elseif outstyle == 3
            tend = tstart - 1.0  # single-step mode
        end

        info::Int = 
        claw1(prob_data,meqn,mwaves,maux,mbc,mx,
              q,aux,xlower,dx,tstart,tend,dtv,cflv,nv,method,mthlim,
              mthbc,
              #work, mwork, # ???
              use_fwaves,
              bc1,rp1,src1,b4step1)

        # check to see if an error occured:
        if info != 0
            println("*** ERROR in claw1 ***  info = $info")
            if info == 1
                println("***   either mx > maxmx or mbc < 2")
            elseif info == 2
                println("***   dt does not divide (tend - tstart)")
                println("***   and dt is fixed since method(1)=0")
            elseif info == 3
                println("***   method(1)=1 and cflv(2) > cflv(1)")
            elseif info == 4
                println("***   mwork is too small")
            elseif info == 11
                println("***   Too many times steps, n > nv(1)")
            elseif info == 12
                println("***   The Courant number is greater than cflv(1)")
                println("***   and dt is fixed since method(1)=0")
            end
 
            return

        end


        dtv[1] = dtv[5]  # use final dt as starting value on next call

        # output solution at this time
        # ------------------------------

        # if outstyle=1 or 2, then nstepout=1 and we output every time
        # we reach this point, since claw1 was called for the entire time
        # increment between outputs.

        # if outstyle=3 then we only output if we have taken nstepout
        # time steps since the last output.

        # iframe is the frame number used to form file names in out1
        iframe::Int = convert(Int, n/nstepout)
        if iframe*nstepout == n
            out1(meqn,mbc,mx,xlower,dx,q,tend,iframe,
                 aux,maux,outaux_always)

###             write(6,601) iframe,tend
###             write(10,1010) tend,info,dtv(3),dtv(4),dtv(5),
###      &           cflv(3),cflv(4),nv(2)
        end
### #
### #        # formats for writing out information about this call to claw:
### #
###   601    format('CLAW1EZ: Frame ',i4,
###      &           ' output files done at time t =',
###      &           d12.4,/)
### #
###  1010    format('tend =',d15.4,/,
###      &       'info =',i5,/,'smallest dt =',d15.4,/,'largest dt =',
###      &       d15.4,/,'last dt =',d15.4,/,'largest cfl =',
###      &         d15.4,/,'last cfl =',d15.4,/,'steps taken =',i4,/)
### #

    end


###   900 continue
###       if (allocated(q))        deallocate(q)
###       if (allocated(aux))      deallocate(aux)
###       if (allocated(work))     deallocate(work)
###       if (allocated(mthlim))   deallocate(mthlim)
###       if (allocated(tout))     deallocate(tout)
###       if (allocated(iout_q))   deallocate(iout_q)
###       if (allocated(iout_aux)) deallocate(iout_aux)

end # claw1ez
