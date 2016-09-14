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

function method(dt::Float64,fc::Int,lc::Int,fm::Int,lm::Int,fs::Int,ls::Int,ifirst::Int,ilast::Int,
                conserved::OffsetArray{Float64,1},
                flux::OffsetArray{Float64,1})
#   ******************************************************************
#   multiply fluxes by dt (ie, compute time integral over cell side)
#   ******************************************************************
    vdt=velocity*dt
    @unsafe for ie=ifirst:ilast+1
        flux[ie]=vdt*conserved[ie-1]
    end
end

