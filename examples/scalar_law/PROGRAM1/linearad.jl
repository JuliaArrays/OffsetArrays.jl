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

function fluxderv(fc::Int,lc::Int,ifirst::Int,ilast::Int,
                  conservd::OffsetArray{Float64,1}, dfdu::OffsetArray{Float64,1})
#   ******************************************************************
#   compute derivative of flux with respect to conserved variable
#   ******************************************************************
    @unsafe for ic=ifirst:ilast
        dfdu[ic]=velocity
    end
end

function riemann(fs::Int,ls::Int,ifirst::Int,ilast::Int,
                 left::OffsetArray{Float64,1},right::OffsetArray{Float64,1},
                 flux::OffsetArray{Float64,1})
#   ******************************************************************
#   the following riemann solver is hard-wired for linear advection
#   with velocity > 0
#   ******************************************************************
    @unsafe for ie=ifirst:ilast+1
        flux[ie]=velocity*left[ie]
    end
end 


