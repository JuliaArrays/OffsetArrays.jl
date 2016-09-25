using OffsetArrays

function src1(meqn::Int, mbc::Int, mx::Int, xlower::Float64, dx::Float64,
              q::OffsetArray{Float64,2},
              maux::Int,
              aux::OffsetArray{Float64,2},
              t::Float64, dt::Float64)
#
#    ! Called to update q by solving source term equation 
#    ! $q_t = \psi(q)$ over time dt starting at time t.
#    !
#    ! This default version does nothing. 
# 
#    implicit none
#    integer, intent(in) :: mbc,mx,meqn,maux
#    real(kind=8), intent(in) :: xlower,dx,t,dt
#    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
#    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)

end # subroutine src1
