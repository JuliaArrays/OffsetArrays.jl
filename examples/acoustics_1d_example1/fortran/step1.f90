

! ===================================================================
    subroutine step1(num_eqn,num_waves,num_ghost,num_aux,mx,q,aux,dx, &
    dt,method,mthlim,cfl,f,wave,s,amdq,apdq,dtdx,use_fwave,rp1)
! ===================================================================

!     # Take one time step, updating q.

!     method(1) = 1   ==>  Godunov method
!     method(1) = 2   ==>  Slope limiter method
!     mthlim(p)  controls what limiter is used in the pth family


!     amdq, apdq, wave, s, and f are used locally:

!     amdq(num_eqn, 1-num_ghost:mx+num_ghost) = left-going flux-differences
!     apdq(num_eqn, 1-num_ghost:mx+num_ghost) = right-going flux-differences
!        e.g. amdq(m,i) = m'th component of A^- \Delta q from i'th Riemann
!                         problem (between cells i-1 and i).

!     wave(num_eqn, num_waves, 1-num_ghost:mx+num_ghost) = waves from solution of
!                                           Riemann problems,
!            wave(m,mw,i) = mth component of jump in q across
!                           wave in family mw in Riemann problem between
!                           states i-1 and i.

!     s(num_waves, 1-num_ghost:mx+num_ghost) = wave speeds,
!            s(m,iw) = speed of wave in family mw in Riemann problem between
!                      states i-1 and i.

!     f(num_eqn, 1-num_ghost:mx+num_ghost) = correction fluxes for second order method
!            f(m,i) = mth component of flux at left edge of ith cell
!     --------------------------------------------------------------------

    implicit double precision (a-h,o-z)
    double precision :: q(num_eqn,1-num_ghost:mx+num_ghost)
    double precision :: aux(num_aux,1-num_ghost:mx+num_ghost)
    double precision :: f(num_eqn,1-num_ghost:mx+num_ghost)
    double precision :: s(num_waves,1-num_ghost:mx+num_ghost)
    double precision :: wave(num_eqn, num_waves,1-num_ghost:mx+num_ghost)
    double precision :: amdq(num_eqn,1-num_ghost:mx+num_ghost)
    double precision :: apdq(num_eqn,1-num_ghost:mx+num_ghost)
    double precision :: dtdx(1-num_ghost:mx+num_ghost)
    integer ::          method(7),mthlim(num_waves)
    logical ::          use_fwave
    logical :: limit
    external :: rp1

!f2py intent(in,out) q
!f2py intent(out) cfl
!f2py intent(in) num_eqn
!f2py intent(in) num_ghost
!f2py intent(in) mx
!f2py optional f, amdq, apdq, dtdx, s, wave


!     # check if any limiters are used:
    limit = .false.
    do 5 mw=1,num_waves
        if (mthlim(mw) > 0) limit = .TRUE. 
    5 END DO

    index_capa = method(6)
    do 10 i=1-num_ghost,mx+num_ghost
        if (index_capa > 0) then
            if (aux(index_capa,i) <= 0.d0) then
                write(6,*) 'Error -- capa must be positive'
                stop
            endif
            dtdx(i) = dt / (dx*aux(index_capa,i))
        else
            dtdx(i) = dt/dx
        endif
    10 END DO



!     # solve Riemann problem at each interface
!     -----------------------------------------

    call rp1(mx,num_eqn,num_waves,num_aux,num_ghost,mx,q,q,aux,aux,wave, &
    s,amdq,apdq)

!     # Modify q for Godunov update:
!     # Note this may not correspond to a conservative flux-differencing
!     # for equations not in conservation form.  It is conservative if
!     # amdq + apdq = f(q(i)) - f(q(i-1)).




          
!     print *,"q before, from fortran",q(24,1)

    do i = 1, mx+1
        ! q(:,i-1) is still in cache from last cycle of i loop, so
        ! update it first
        do m = 1, num_eqn
            q(m,i-1) = q(m,i-1) - dtdx(i-1)*amdq(m,i)
        end do
        do m = 1, num_eqn
            q(m,i) = q(m,i) - dtdx(i)*apdq(m,i)
        end do
    end do


!     # compute maximum wave speed:
    cfl = 0.d0
    do i=1,mx+1
        do mw=1,num_waves
        !          # if s>0 use dtdx(i) to compute CFL,
        !          # if s<0 use dtdx(i-1) to compute CFL:
            cfl = dmax1(cfl, dtdx(i)*s(mw,i), -dtdx(i-1)*s(mw,i))
        end do
    end do

    if (method(2) == 1) go to 900

!     # compute correction fluxes for second order q_{xx} terms:
!     ----------------------------------------------------------

!      # apply limiter to waves:
    if (limit) call limiter(mx,num_eqn,num_waves,num_ghost,mx, &
    wave,s,mthlim)

    if (use_fwave.eqv. .FALSE. ) then
        do i=1,mx+1
            do m = 1,num_eqn
                f(m,i) = 0.d0
            end do
            dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
            do mw=1,num_waves
                do m=1,num_eqn
                    f(m,i) = f(m,i) + 0.5d0 * dabs(s(mw,i)) &
                    * (1.d0 - dabs(s(mw,i))*dtdxave) * wave(m,mw,i)
                end do
            end do
        end do
    else
        do i=1,mx+1
            do m = 1,num_eqn
                f(m,i) = 0.d0
            end do
            dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
            do mw=1,num_waves
                do m=1,num_eqn
                    f(m,i) = f(m,i) + 0.5d0 * dsign(1.d0,s(mw,i)) &
                    * (1.d0 - dabs(s(mw,i))*dtdxave) * wave(m,mw,i)
                end do
            end do
        end do
    endif



!     # update q by differencing correction fluxes
!     ============================================

!     # (Note:  Godunov update has already been performed above)

    do i = 1, mx+1
        do m = 1, num_eqn
            q(m,i) = q(m,i) - dtdx(i) * (f(m,i+1) - f(m,i))
        end do
    end do

    900 continue
    return
    end subroutine step1
