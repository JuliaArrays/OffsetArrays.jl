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

function consdiff(fic::Int,lac::Int,fif::Int,laf::Int,fix::Int,lax::Int,
                  ifirst::Int,ilast::Int,
                  x::OffsetArray{Float64,1},
                  flux::OffsetArray{Float64,1},
                  conservd::OffsetArray{Float64,1})
#   ******************************************************************
#   update conservd to new time
#   ******************************************************************
    @unsafe for ic=ifirst:ilast
        conservd[ic] -= (flux[ic+1]-flux[ic]) / (x[ic+1]-x[ic])
    end
end

function stabledt(fc::Int,lc::Int,fm::Int,lm::Int,ifirst::Int,ilast::Int,
                  conserved::OffsetArray{Float64,1},
                  lambda_cell::OffsetArray{Float64,1},
                  x::OffsetArray{Float64,1})
#   ******************************************************************
#   compute stable timestep
#   ******************************************************************
    sdt=huge
    @unsafe for ic=ifirst:ilast
        abs_lambda=abs(lambda_cell[ic])
        dx = x[ic+1]-x[ic]
        if abs_lambda > roundoff*dx
            sdt=min(sdt,dx/abs_lambda)
        end
    end
    @assert sdt > 0.0
    sdt
end


