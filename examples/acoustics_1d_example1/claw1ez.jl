
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

    unit10 = open("fort.info", "w")
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

    outstyle = readint(r)
 
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

    elseif outstyle == 3
        nstepout  = readint(r)
        nstop     = readint(r)
        output_t0 = readbool(r)       # Not currently used
        nout      = nstop

    else
        println("nout: $nout tfinal: $tfinal output_t0: $output_t0")
        println("*** Unrecognized output style $outstyle")
        error("*** Exiting claw1ez")
    end

    output_format = readint(r)    # Not used yet

    iout_q = readfltarray(r, meqn)
    println("iout_q: $iout_q")

#        iout_aux

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

    dt_variable::Bool = readbool(r)
    method[1] = dt_variable ? 1 : 0

    method[2] = readint(r) # Order

    # method(3) (transverse order) not used in 1D
    # No dimensional splitting in 1D

    method[4] = readint(r) # Verbosity
    method[5] = readint(r) # Source term splitting style
    method[6] = readint(r) # Index into aux of capacity function
    method[7] = maux       # Number of aux variables

    println("method: $method")

    use_fwaves::Bool = readbool(r)

    mthlim::Array{Int, 1} = readintarray(r, mwaves)
    println("mthlim: $mthlim")

    mbc      = readint(r)
    mthbc[1] = readint(r)
    mthbc[2] = readint(r)

    println("mbc: $mbc mthbc: $mthbc")


    # No restart in 1D, and there's nothing after the restart info, so
    # just close the file now.
    close_reader(r) 

    # Check consistency for periodic BCs
    if (mthbc[1] == 2 && mthbc[2] != 2) ||
       (mthbc[2] == 2 && mthbc[1] != 2)
          error("*** ERROR ***  periodic boundary conditions require mthbc(1) and mthbc(2) BOTH be set to 2")
    end

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

    q = OffsetArray(Float64, 1:meqn, 1-mbc:mx+mbc)

    # set initial conditions:

    qinit(prob_data, meqn, mbc, mx, xlower, dx, q, maux, aux)

    # output initial data
    out1(meqn,mbc,mx,xlower,dx,q,t0,0,aux,maux,
         (outaux_init_only || outaux_always) )
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

            @printf("CLAW1EZ: Frame %4d output files done at time t = %12.4e\n", iframe,tend)
            @printf(unit10,"tend = %15.4e \ninfo =%5d\nsmallest dt =%15.4e\nlargest dt =%15.4e\n",
                    tend,info,dtv[3],dtv[4])
            @printf(unit10,"last dt =%15.4e\nlargest cfl =%15.4e\nlast cfl =%15.4e\nsteps taken =%4d\n",
                    dtv[5],cflv[3],cflv[4],nv[2])

        end
    end
    Base.close(unit10)

end # claw1ez
