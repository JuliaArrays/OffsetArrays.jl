using OffsetArrays

#
#
#=========================================================================#
function bc1(meqn::Int, mbc::Int, mx::Int, xlower::Float64, dx::Float64,
             q::OffsetArray{Float64,2}, maux::Int, aux::OffsetArray{Float64,2},
             t::Float64, dt::Float64, mthbc::Array{Int,1})
#=========================================================================#
#
#     # Standard boundary condition choices for claw2
#
#     # At each boundary  k = 1 (left),  2 (right):
#     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
#     #            =  1  for zero-order extrapolation
#     #            =  2  for periodic boundary coniditions
#     #            =  3  for solid walls, assuming this can be implemented
#     #                  by reflecting the data about the boundary and then
#     #                  negating the 2'nd component of q.
#     ------------------------------------------------
#
#     # Extend the data from the computational region
#     #      i = 1, 2, ..., mx2
#     # to the virtual cells outside the region, with
#     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
#
##      implicit double precision (a-h,o-z)
##      dimension q(meqn,1-mbc:mx+mbc)
##      dimension aux(maux,1-mbc:mx+mbc)
##
##      dimension mthbc(2)

#
#
#-------------------------------------------------------
#     # left boundary:
#-------------------------------------------------------
#      go to (100,110,120,130) mthbc(1)+1
#
    if mthbc[1] == 0
#  100 continue

        # user-specified boundary conditions go here in place of error output
        println("*** ERROR *** mthbc(1)=0 and no BCs specified in bc1")

    elseif mthbc[1] == 1

#  110 continue
        # zero-order extrapolation:
        for ibc=1:mbc
            for m=1:meqn
                q[m,1-ibc] = q[m,1]
            end
        end
#      go to 199

    elseif mthbc[1] == 2 
#  120 continue
        # periodic:
        for ibc=1:mbc
            for m=1:meqn
                q[m,1-ibc] = q[m,mx+1-ibc]
            end
        end
#      go to 199

    elseif mthbc[1] == 3 
#  130 continue
#     # solid wall (assumes 2'nd component is velocity or momentum in x):
      for ibc=1:mbc
         for m=1:meqn
            q[m,1-ibc] = q[m,ibc]
         end
         q[2,1-ibc] = -q[2,ibc]
      end
#     # negate the normal velocity:
#      go to 199

#  199 continue
    end # if mthbc[1] == ...

#
#-------------------------------------------------------
#     # right boundary:
#-------------------------------------------------------
#      go to (200,210,220,230) mthbc(2)+1
#
    if mthbc[2] == 0
#  200 continue
        # user-specified boundary conditions go here in place of error output
        println("*** ERROR *** mthbc(2)=0 and no BCs specified in bc2")
        # go to 299

    elseif mthbc[2] == 1
#  210 continue
        # zero-order extrapolation:
        for ibc=1:mbc
            for m=1:meqn
                q[m,mx+ibc] = q[m,mx]
            end
        end
#      go to 299

    elseif mthbc[2] == 2 
#  220 continue
#     # periodic:  
        for ibc=1:mbc
            for m=1:meqn
                q[m,mx+ibc] = q[m,ibc]
            end
        end
#      go to 299

    elseif mthbc[2] == 3 
#  230 continue
#     # solid wall (assumes 2'nd component is velocity or momentum in x):
      for ibc=1:mbc
         for m=1:meqn
            q[m,mx+ibc] = q[m,mx+1-ibc]
         end
         q[2,mx+ibc] = -q[2,mx+1-ibc]
      end
#      go to 299

    end # if mthbc[2] == ...
#  299 continue
#

end

