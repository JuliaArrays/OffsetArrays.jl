function b4step1(mbc::Int, mx::Int, meqn::Int, q, xloweri::Float64, dxi::Float64, t::Float64, dt::Float64, mauxi::Int, aux)
#
#    ! Called before each call to step1.
#    ! Use to set time-dependent aux arrays or perform other tasks.
#    !
#    ! This default version does nothing. 
# 
#    implicit none
#    integer, intent(in) :: mbc,mx,meqn,maux
#    real(kind=8), intent(in) :: xlower,dx,t,dt
#    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)
#    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc)
#
end # subroutine b4step1

