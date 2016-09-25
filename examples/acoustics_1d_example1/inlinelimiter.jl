#
# Inlined version of limiter code
#
#     =====================================================
#     # Apply a limiter to the waves.
#     # The limiter is computed by comparing the 2-norm of each wave with
#     # the projection of the wave from the interface to the left or
#     # right onto the current wave.  For a linear system this would
#     # correspond to comparing the norms of the two waves.  For a
#     # nonlinear problem the eigenvectors are not colinear and so the
#     # projection is needed to provide more limiting in the case where the
#     # neighboring wave has large norm but points in a different direction
#     # in phase space.
#
#     # The specific limiter used in each family is determined by the
#     # value of the corresponding element of the array mthlim, as used in
#     # the function philim.
#     # Note that a different limiter may be used in each wave family.
#
#     # dotl and dotr denote the inner product of wave with the wave to
#     # the left or right.  The norm of the projections onto the wave are then
#     # given by dotl/wnorm2 and dotr/wnorm2, where wnorm2 is the 2-norm
#     # of wave.

using OffsetArrays

function limiter(maxm::Int, num_eqn::Int, num_waves::Int, num_ghost::Int, mx::Int, wave::OffsetArray{Float64}, s::OffsetArray{Float64}, mthlim::Array{Int,1})
#
#    # Arguments
#    integer, intent(in) :: maxm, num_eqn, num_waves, num_ghost, mx
#    real(kind=8), intent(in out) :: wave(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
#    real(kind=8), intent(in) :: s(num_waves, 1-num_ghost:maxm+num_ghost)
#    integer, intent(in) :: mthlim(num_waves)
#
#    ! Local storage
#    integer :: m, mw, i
#    real(kind=8) :: r, c, wlimiter, wnorm2, dotl
#    real(kind=8), dimension(num_waves) :: dotr
    dotr = Array(Float64, num_waves)

    dotr[:] = 0.0

    # x_loop:
    for i = 0:mx+1

        # wave_loop:
        for mw = 1:num_waves
            if mthlim[mw] == 0
                continue # wave_loop
            end

            # Construct dot products
            wnorm2 = 0.0
            dotl = dotr[mw]
            dotr[mw] = 0.0
            for m=1:num_eqn
                wnorm2 = wnorm2 + wave[m,mw,i]^2
                dotr[mw] = dotr[mw] + wave[m,mw,i]*wave[m,mw,i+1]
            end

            # Skip this loop if it's on the boundary or the size of the wave is
            # zero (but still want dot products to be initialized above)
            if i == 0
                continue # cycle wave_loop
            end
            if wnorm2 == 0.0
                continue # cycle wave_loop
            end

            # Compute ratio of this wave's strength to upwind wave's strength
            if s[mw,i] > 0.0
                r = dotl / wnorm2
            else
                r = dotr[mw] / wnorm2
            end

            # Compute value of limiter function
            # Minmod
            if     mthlim[mw] == 1
                wlimiter = max(0.0, min(1.0, r))

            # Superbee
            elseif mthlim[mw] == 2
                wlimiter = max(0.0, min(1.0, 2.0*r), min(2.0, r))

                # Van Leer
            elseif mthlim[mw] == 3 
                wlimiter = (r + abs(r)) / (1.0 + abs(r))

                # Monotonized - Centered
            elseif mthlim[mw] == 4 
                c = (1.0 + r)/2.0
                wlimiter = max(0.0, min(c, 2.0, 2.0*r))

                # Beam Warming
            elseif mthlim[mw] == 5
                wlimiter = r

            else
                error("Invalid limiter method.")

            end # if

            # Apply resulting limit
            # broken below => expand into loop
            #wave[..,mw,i] = wlimiter * wave[..,mw,i]
            for k = 1:num_eqn
                wave[k,mw,i] = wlimiter * wave[k,mw,i]
            end
 
        end # do wave_loop
    end # do x_loop

end # subroutine limiter
