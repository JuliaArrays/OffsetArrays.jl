using OffsetArrays
using SetProb

# ============================================================================
function rp1(parms::CParam,
             maxm::Int, meqn::Int, mwaves::Int, maux::Int, mbc::Int, mx::Int,
             ql::OffsetArray{Float64,2}, qr::OffsetArray{Float64,2},
             auxl::OffsetArray{Float64}, auxr::OffsetArray{Float64},
             wave::OffsetArray{Float64,3}, s::OffsetArray{Float64,2},
             amdq::OffsetArray{Float64,2}, apdq::OffsetArray{Float64,2})
# ============================================================================

# subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)

# Riemann solver for the acoustics equations in 1d,

# On input, ql contains the state vector at the left edge of each cell
#           qr contains the state vector at the right edge of each cell

# On output, wave contains the waves,
#            s the speeds,
#
#            amdq = A^- Delta q,
#            apdq = A^+ Delta q,
#                   the decomposition of the flux difference
#                       f(qr(i-1)) - f(ql(i))
#                   into leftgoing and rightgoing parts respectively.
#

# Note that the i'th Riemann problem has left state qr(i-1,:)
#                                    and right state ql(i,:)
# From the basic clawpack routines, this routine is called with ql = qr


#    dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
#    dimension    s(mwaves,1-mbc:maxm+mbc)
#    dimension   ql(meqn, 1-mbc:maxm+mbc)
#    dimension   qr(meqn, 1-mbc:maxm+mbc)
#    dimension apdq(meqn, 1-mbc:maxm+mbc)
#    dimension amdq(meqn, 1-mbc:maxm+mbc)

#     local arrays
#     ------------
#    dimension delta(2)
    delta = Array(Float64, 2)

# density, bulk modulus, and sound speed, and impedence of medium:
# (should be set in setprob.f)
#    common /cparam/ rho,bulk,cc,zz
 
    # split the jump in q at each interface into waves
 
    # find a1 and a2, the coefficients of the 2 eigenvectors:

    for i = 2-mbc : mx+mbc
        delta[1] = ql[1,i] - qr[1,i-1]
        delta[2] = ql[2,i] - qr[2,i-1]
        a1 = (-delta[1] + parms.zz*delta[2]) / (2.0*parms.zz)
        a2 = ( delta[1] + parms.zz*delta[2]) / (2.0*parms.zz)
     
        # Compute the waves.
 
        wave[1,1,i] = -a1*parms.zz
        wave[2,1,i] = a1
        s[1,i] = -parms.cc

        wave[1,2,i] = a2*parms.zz
        wave[2,2,i] = a2
        s[2,i] = parms.cc

    end

    # compute the leftgoing and rightgoing flux differences:
    # Note s(1,i) < 0   and   s(2,i) > 0.

    for m = 1 : meqn
        for i = 2-mbc : mx+mbc
            amdq[m,i] = s[1,i]*wave[m,1,i]
            apdq[m,i] = s[2,i]*wave[m,2,i]
        end 
    end

end #subroutine rp1
