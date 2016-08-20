using BenchmarkTools
using OffsetArrays

const dim = 1000

x = Array(Float64, 2*dim)
y = OffsetArray(Float64, -dim + 1 : dim)

fillx(x) = for i in indices(x,1); x[i] = i; end
updatex(x) = for i in indices(x,1); x[i] = x[i] + i; end
updatex_eachindex(x) = for i in eachindex(x); x[i] = x[i] + i; end

@show @benchmark fillx(x)
@show @benchmark fillx(y)
@show @benchmark updatex(x)
@show @benchmark updatex(y)
@show @benchmark updatex_eachindex(x)
@show @benchmark updatex_eachindex(y)
