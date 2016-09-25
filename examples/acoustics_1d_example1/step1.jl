using OffsetArrays

include("rp1_acoustics.jl")
include("inlinelimiter.jl")

# ===================================================================
function step1(parms::CParam,num_eqn::Int, num_waves::Int, num_ghost::Int, num_aux::Int, mx::Int,
               q::OffsetArray{Float64}, aux::OffsetArray{Float64},
               dx::Float64, dt::Float64, method::Array{Int,1}, mthlim::Array{Int},
               f::OffsetArray{Float64}, wave::OffsetArray{Float64}, s::OffsetArray{Float64},
               amdq::OffsetArray{Float64}, apdq::OffsetArray{Float64}, dtdx::OffsetArray{Float64},
               use_fwave::Bool, rp1)
# ===================================================================

    # Take one time step, updating q.

    # method(1) = 1   ==>  Godunov method
    # method(1) = 2   ==>  Slope limiter method
    # mthlim(p)  controls what limiter is used in the pth family
#
#
#     amdq, apdq, wave, s, and f are used locally:
#
#     amdq(num_eqn, 1-num_ghost:mx+num_ghost) = left-going flux-differences
#     apdq(num_eqn, 1-num_ghost:mx+num_ghost) = right-going flux-differences
#        e.g. amdq(m,i) = m'th component of A^- \Delta q from i'th Riemann
#                         problem (between cells i-1 and i).
#
#     wave(num_eqn, num_waves, 1-num_ghost:mx+num_ghost) = waves from solution of
#                                           Riemann problems,
#            wave(m,mw,i) = mth component of jump in q across
#                           wave in family mw in Riemann problem between
#                           states i-1 and i.
#
#     s(num_waves, 1-num_ghost:mx+num_ghost) = wave speeds,
#            s(m,iw) = speed of wave in family mw in Riemann problem between
#                      states i-1 and i.
#
#     f(num_eqn, 1-num_ghost:mx+num_ghost) = correction fluxes for second order method
#            f(m,i) = mth component of flux at left edge of ith cell
#     --------------------------------------------------------------------
#
#    implicit double precision (a-h,o-z)
#    double precision :: q(num_eqn,1-num_ghost:mx+num_ghost)
#    double precision :: aux(num_aux,1-num_ghost:mx+num_ghost)
#    double precision :: f(num_eqn,1-num_ghost:mx+num_ghost)
#    double precision :: s(num_waves,1-num_ghost:mx+num_ghost)
#    double precision :: wave(num_eqn, num_waves,1-num_ghost:mx+num_ghost)
#    double precision :: amdq(num_eqn,1-num_ghost:mx+num_ghost)
#    double precision :: apdq(num_eqn,1-num_ghost:mx+num_ghost)
#    double precision :: dtdx(1-num_ghost:mx+num_ghost)


    cfl::Float64 = 0.0

#    dtdx = OffsetArray(Float64, 1-num_ghost:mx+num_ghost)
#    integer ::          method(7),mthlim(num_waves)
#    logical ::          use_fwave
#    logical :: limit
#    external :: rp1

#f2py intent(in,out) q
#f2py intent(out) cfl
#f2py intent(in) num_eqn
#f2py intent(in) num_ghost
#f2py intent(in) mx
#f2py optional f, amdq, apdq, dtdx, s, wave


    # check if any limiters are used:
    limit = false
    for mw=1:num_waves
        if mthlim[mw] > 0
            limit = true
            break
        end
    end

    index_capa = method[6]
    for i=1-num_ghost:mx+num_ghost
        if index_capa > 0
            if aux[index_capa,i] <= 0.0
                error("capa must be positive")
            end
            dtdx[i] = dt / (dx*aux[index_capa,i])
        else
            dtdx[i] = dt / dx
        end
    end


    # solve Riemann problem at each interface
    # -----------------------------------------

    rp1(parms,mx,num_eqn,num_waves,num_aux,num_ghost,mx,q,q,aux,aux,wave,s,amdq,apdq)


    # Modify q for Godunov update:
    # Note this may not correspond to a conservative flux-differencing
    # for equations not in conservation form.  It is conservative if
    # amdq + apdq = f(q(i)) - f(q(i-1)).
    
    # print *,"q before, from fortran",q(24,1)

    for i = 1 : mx+1
        # q(:,i-1) is still in cache from last cycle of i loop, so
        # update it first
        for m = 1 : num_eqn
            q[m,i-1] = q[m,i-1] - dtdx[i-1]*amdq[m,i]
        end
        for m = 1 : num_eqn
            q[m,i] = q[m,i] - dtdx[i]*apdq[m,i]
        end
    end

    # compute maximum wave speed:
    cfl = 0.0
    for i=1:mx+1
        for mw=1:num_waves
            # if s>0 use dtdx(i) to compute CFL,
            # if s<0 use dtdx(i-1) to compute CFL:
            #cfl = dmax1(cfl, dtdx(i)*s(mw,i), -dtdx(i-1)*s(mw,i))
            #cfl = max(cfl, dtdx[i]*s[mw,i], -dtdx[i-1]*s[mw,i])
            s1 = s[mw,i]
            cfl = s1 > 0.0 ? s1*dtdx[i] : -s1*dtdx[i-1]
        end
    end

    if method[2] == 1
        return
    end

    # compute correction fluxes for second order q_{xx} terms:
    #     ----------------------------------------------------------

    # apply limiter to waves:
    if limit
        limiter(mx,num_eqn,num_waves,num_ghost,mx,wave,s,mthlim)
    end

    if !use_fwave
        for i=1:mx+1
            for m = 1:num_eqn
                f[m,i] = 0.0
            end
            dtdxave = 0.5 * (dtdx[i-1] + dtdx[i])
            for mw=1:num_waves
                for m=1:num_eqn
                    f[m,i] += 0.5 * abs(s[mw,i]) * (1.0 - abs(s[mw,i])*dtdxave) * wave[m,mw,i]
                end
            end
        end
    else
        for i=1:mx+1
            for m = 1:num_eqn
                f[m,i] = 0.0
            end
            dtdxave = 0.5 * (dtdx[i-1] + dtdx[i])
            for mw=1:num_waves
                for m=1:num_eqn
                    f[m,i] += 0.5 * sign(s[mw,i]) * (1.0 - abs(s[mw,i])*dtdxave) * wave[m,mw,i]
                end
            end
        end
    end



    # update q by differencing correction fluxes
    # ============================================

    # (Note:  Godunov update has already been performed above)

    for i = 1: mx+1
        for m = 1: num_eqn
            q[m,i] -= dtdx[i] * (f[m,i+1] - f[m,i])
        end
    end

    return cfl

end # subroutine step1
