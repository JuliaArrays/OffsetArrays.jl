using BenchmarkTools
using OffsetArrays

const dim = 1000

x = Array(Float64, 2*dim)
y = OffsetArray(Float64, -dim + 1 : dim)

fillx() = for i = 1:2*dim x[i] = i end
filly() = for i = -dim+1:dim y[i] = i end
updatex() = for i = 1:2*dim x[i] = i end
updatey() = for i = -dim+1:dim y[i] = i end

@benchmark fillx() 
@benchmark filly()
@benchmark updatex() 
@benchmark updatey()

