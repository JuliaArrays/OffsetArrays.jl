#***********************************************************************
#  Copyright 2006 John A. Trangenstein
#
#  This software is made available for research and instructional use 
#  only. 
#  You may copy and use this software without charge for these 
#  non-commercial purposes, provided that the copyright notice and 
#  associated text is reproduced on all copies.  
#  For all other uses (including distribution of modified versions), 
#  please contact the author at
#    John A. Trangenstein
#    Department of Mathematics
#    Duke University
#    Durham, NC 27708-0320
#    USA
#  or
#    johnt@math.duke.edu
#  
#  This software is made available "as is" without any assurance that it
#  is completely correct, or that it will work for your purposes.  
#  Use the software at your own risk.
#***********************************************************************

using OffsetArrays

const roundoff  =  1.0-14
const small     =  1.0-20
const huge      =  1.0e300
const undefind  =  1.0e300
const lnrndoff  =  14.0e0
const lnsmall   =  20.0e0

# problem-specific parameters:
const jump      =  0.0
const x_left    = -0.2
const x_right   =  1.0
const statelft  =  2.0
const statergt  =  0.0
const velocity  =  1.0

include("riemprob.jl")
include("linearad.jl")
include("consdiff.jl")
include("upwind.jl")

function main()
    const ncells = 10000

    u    = OffsetArray(Float64, -2:ncells+1)
    x    = OffsetArray(Float64,  0:ncells)
    flux = OffsetArray(Float64,  0:ncells)
    dfdu = OffsetArray(Float64, -2:ncells+1)

    const nsteps   =  100
    #const nsteps   =  10000
    const tmax     =  0.8
    const cfl      =  0.9

    # array bounds:
    const fc=-2
    const lc=ncells+1
    const fm=0
    const lm=ncells-1
    const fs=0
    const ls=ncells-1
    const ifirst=0
    const ilast=ncells-1

    initsl(ncells,fc,lc,fm,lm,ifirst,ilast, u,x)

    #println("x : $x u : $u")
    bcmesh(fm,lm,ncells, x)
    bccells(fc,lc,ncells, u)
    fluxderv(fc,lc,fc,lc, u, dfdu)

    istep=0
    t=0.0
    while istep < nsteps && t < tmax
        bccells(fc,lc,ncells, u)
        fluxderv(fc,lc,fc,lc, u, dfdu)
        dt=cfl*stabledt(fc,lc,fm,lm,ifirst,ilast, u, dfdu,x)
        method(dt,fc,lc,fm,lm,fs,ls,ifirst,ilast, u, flux)
        consdiff(fc,lc,fs,ls,fm,lm,ifirst,ilast, x,flux, u)

        t=t+dt
        istep=istep+1
    end

    # write final results (plot later)
    #@unsafe for ic=0:ncells-1
    for ic=0:ncells-1
        xc = (x[ic]+x[ic+1])*0.5
        uc = u[ic]
        @printf("%e %e\n",xc,uc)
    end


end # main

main()


