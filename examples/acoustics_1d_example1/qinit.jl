using OffsetArrays
using SetProb

function qinit(parm::CParam, meqn::Int, mbc::Int, mx::Int, xlower::Float64, dx::Float64, q::OffsetArray{Float64}, maux::Int, aux::OffsetArray{Float64})

    # Set initial conditions for the q array.
    # This default version prints an error message since it should
    # not be used directly.  Copy this to an application directory and
    # loop over all grid cells to set values of q(1:meqn, 1:mx).

    
    # integer, intent(in) :: meqn,mbc,mx,maux
    # real(kind=8), intent(in) :: xlower,dx
    # real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    # real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    # integer :: i
    # real(kind=8) :: beta, xcell
    # common /cqinit/ beta
 
 
    for i = 1 : mx
        xcell = xlower + (i-0.5)*dx
        q[1,i] = exp(-parm.beta * (xcell-0.3)^2)  
        q[2,i] = 0.0
        #println("q[1,$i] = $(q[1,i])")
        #println("q[2,$i] = $(q[2,i])")
    end

end # qinit

