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
function main()

    const ncells=10000

    #   integer nsteps
    #   double precision cfl,tmax
    #   double precision jump,x_left,x_right,statelft,statergt,velocity
    #   double precision
    #  &  u(-2:ncells+1),
    #  &  x(0:ncells),
    #  &  flux(0:ncells)

    u    = Array(Float64, ncells+4)
    x    = Array(Float64, ncells+1)
    flux = Array(Float64, ncells+1)

   #   integer fc,lc,ifirst,ilast
   #   integer ic,ie,ijump,istep
   #   double precision dt,dx,frac,mindx,t,vdt
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #     problem-specific parameters:
    tic()
    const jump     =  0.0
    const x_left   = -0.2
    const x_right  =  1.0
    const statelft =  2.0
    const statergt =  0.0
    const velocity =  1.0

    const nsteps   =  10000
    const tmax     =  0.8
    const cfl      =  0.9

    # array bounds:
    const fc=-2
    const lc=ncells+1
    const ifirst=0
    const ilast=ncells-1

    # uniform mesh:
    dx=(x_right-x_left)/ncells
    for ie in ifirst:ilast+1
        x[ie+1]=x_left+ie*dx
    end

    # initial values for diffential equation:
    ijump=max(ifirst-1,min(convert(Int,round(ncells*(jump-x_left)/(x_right-x_left))),ilast+1))

    # left state to left of jump
    for ic=ifirst:ijump-1
        u[ic+3]=statelft
    end

    # volume-weighted average in cell containing jump
    frac=(jump-x_left-ijump*dx)/(x_right-x_left)
    u[ijump+3]=statelft*frac+statergt*(1.0-frac)

    # right state to right of jump
    for ic=ijump+1:ilast
        u[ic+3]=statergt
    end

    # stable timestep (independent of time for linear advection):
    mindx=1.0e300
    for ic=ifirst:ilast
        mindx=min(mindx,x[ic+2]-x[ic+1])
    end
    dt=cfl*mindx/abs(velocity)

    istep=0
    t=0.0

    # loop over timesteps
    while istep < nsteps && t < tmax
        # right boundary condition: outgoing wave
        for ic=ncells:lc
            u[ic+3]=u[ncells+2]
        end

        # left boundary condition: specified value
        for ic=fc:-1
          u[ic+3]=statelft
        end

        # upwind fluxes times dt (ie, flux time integral over cell side)
        # assumes velocity > 0
        vdt=velocity*dt
        for ie=ifirst:ilast+1
          flux[ie+1]=vdt*u[ie+2]
        end

        # conservative difference
        for ic=ifirst:ilast
          u[ic+3] -= (flux[ic+2]-flux[ic+1]) / (x[ic+2]-x[ic+1])
        end

        # update time and step number
        t=t+dt
        istep=istep+1
    end

    # write final results (plot later)
    for ic=0:ncells-1
        xc = (x[ic+1]+x[ic+2])*0.5
        uc = u[ic+3]
        @printf("%e %e\n",xc,uc)
    end
    toc()
end # main

main()

# @profile main()
# Profile.print()

