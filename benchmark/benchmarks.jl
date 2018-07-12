using BenchmarkTools
using OffsetArrays

const dim = 1000

x = Array{Float64}(undef, 2*dim)
y = OffsetArray{Float64}(undef, -dim + 1 : dim)
x2d = Array{Float64}(undef, 2*dim, 2*dim)
y2d = OffsetArray{Float64}(undef, -dim + 1 : dim, -dim + 1 : dim)

fill(x) = for i in axes(x,1); x[i] = i; end
fill2d(x) = for j in axes(x,2); for i in axes(x,1); x[i,j] = i + j; end; end
update(x) = for i in axes(x,1); x[i] = x[i] + i; end
update2d(x) = for j in axes(x,2); for i in axes(x,1); x[i,j] = x[i,j] + i + j; end; end
update_eachindex(x) = for i in eachindex(x); x[i] = x[i] + i; end

unsafe_fill(x) = @inbounds(for i in axes(x,1); x[i] = i; end)
unsafe_fill2d(x) = @inbounds(for j in axes(x,2); for i in axes(x,1); x[i,j] = i + j; end; end)
unsafe_update(x) = @inbounds(for i in axes(x,1); x[i] = x[i] + i; end)
unsafe_update2d(x) = @inbounds(for j in axes(x,2); for i in axes(x,1); x[i,j] = x[i,j] + i + j; end; end)
unsafe_update_eachindex(x) = @inbounds(for i in eachindex(x); x[i] = x[i] + i; end)

@show @benchmark fill(x)
@show @benchmark fill(y)
@show @benchmark fill2d(x2d)
@show @benchmark fill2d(y2d)
@show @benchmark update(x)
@show @benchmark update(y)
@show @benchmark update2d(x2d)
@show @benchmark update2d(y2d)
@show @benchmark update_eachindex(x)
@show @benchmark update_eachindex(y)

@show @benchmark unsafe_fill(x)
@show @benchmark unsafe_fill(y)
@show @benchmark unsafe_fill2d(x2d)
@show @benchmark unsafe_fill2d(y2d)
@show @benchmark unsafe_update(x)
@show @benchmark unsafe_update(y)
@show @benchmark unsafe_update2d(x2d)
@show @benchmark unsafe_update2d(y2d)
@show @benchmark unsafe_update_eachindex(x)
@show @benchmark unsafe_update_eachindex(y)
