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

function initsl(ncells::Int, fc::Int, lc::Int, fx::Int, lx::Int, ifirst::Int, ilast::Int,
                conservd::OffsetArray{Float64,1},
                x::OffsetArray{Float64,1})
#   ******************************************************************
#   interior values only; others defined by boundary conditions
#   ******************************************************************
    dx=(x_right-x_left)/ncells
    @unsafe for ie=ifirst:ilast+1
        x[ie]=x_left+ie*dx
    end
    ijump=max(ifirst-1,min(convert(Int,round(ncells*(jump-x_left)/(x_right-x_left))),ilast+1))
    @unsafe for ic=ifirst:ijump-1
        conservd[ic]=statelft
    end
    frac=(jump-x_left-ijump*dx)/(x_right-x_left)
    conservd[ijump]=statelft*frac+statergt*(1.0-frac)
    @unsafe for ic=ijump+1:ilast
        conservd[ic]=statergt
    end
end

function bcmesh(fim::Int,lam::Int,ncells::Int, x::OffsetArray{Float64,1})
#   ******************************************************************
#   right side 
#   ******************************************************************
    if lam >= ncells
        it=ncells
        dx=x(it)-x(it-1)
        @unsafe for ie=it+1:lam+1
            x[ie]=x[it]+dx*(ie-it)
        end
    end
#   ******************************************************************
#   left side 
#   ******************************************************************
    if fim < 0
        dx=x(1)-x(0)
        @unsafe for ie=fim:-1
            x[ie]=x[0]+dx*(ie)
        end
    end      
end

function bccells(fic::Int,lac::Int,ncells::Int, conservd::OffsetArray{Float64,1})
#   ******************************************************************
#   right side 
#   ******************************************************************
    if lac >= ncells
        # println("treating right side")
        @unsafe for ic=ncells:lac
            # outgoing wave:
            conservd[ic]=conservd[ncells-1]
        end  
    end
#     ******************************************************************
#     left side 
#     ******************************************************************
    if fic < 0
        # println("treating left side")
        @unsafe for ic=fic:-1
            # outgoing wave:
            # conservd[ic]=conservd[0]

            # specified value
            conservd[ic]=statelft
        end  
    end
end

