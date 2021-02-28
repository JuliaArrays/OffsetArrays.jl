using BenchmarkTools
using OffsetArrays

const dim = 1000

x = Array{Float64}(undef, 2*dim);
y = OffsetArray{Float64}(undef, -dim + 1 : dim);
x2d = Array{Float64}(undef, 2*dim, 2*dim);
y2d = OffsetArray{Float64}(undef, -dim + 1 : dim, -dim + 1 : dim);

r = OffsetVector(1:dim, 0);
s = OffsetVector(1:dim, 0);

fill1d(x) = for i in axes(x,1); x[i] = i; end
fill2d(x) = for j in axes(x,2); for i in axes(x,1); x[i,j] = i + j; end; end
update(x) = for i in axes(x,1); x[i] = x[i] + i; end
update2d(x) = for j in axes(x,2); for i in axes(x,1); x[i,j] = x[i,j] + i + j; end; end
update_eachindex(x) = for i in eachindex(x); x[i] = x[i] + i; end

unsafe_fill(x) = @inbounds(for i in axes(x,1); x[i] = i; end)
unsafe_fill2d(x) = @inbounds(for j in axes(x,2); for i in axes(x,1); x[i,j] = i + j; end; end)
unsafe_update(x) = @inbounds(for i in axes(x,1); x[i] = x[i] + i; end)
unsafe_update2d(x) = @inbounds(for j in axes(x,2); for i in axes(x,1); x[i,j] = x[i,j] + i + j; end; end)
unsafe_update_eachindex(x) = @inbounds(for i in eachindex(x); x[i] = x[i] + i; end)

vectorlinearindexing(a, ax) = a[ax]
vectorCartesianindexing(a, ax1, ax2) = a[ax1, ax2]
nestedvectorlinearindexing(a, ax1, ax2) = a[ax1[ax2]]

macro showbenchmark(ex)
	print(ex, " : ")
	quote
		println(@benchmark esc($ex))
	end
end

@showbenchmark fill1d(x)
@showbenchmark fill1d(y)
@showbenchmark fill2d(x2d)
@showbenchmark fill2d(y2d)
@showbenchmark update(x)
@showbenchmark update(y)
@showbenchmark update2d(x2d)
@showbenchmark update2d(y2d)
@showbenchmark update_eachindex(x)
@showbenchmark update_eachindex(y)

@showbenchmark unsafe_fill(x)
@showbenchmark unsafe_fill(y)
@showbenchmark unsafe_fill2d(x2d)
@showbenchmark unsafe_fill2d(y2d)
@showbenchmark unsafe_update(x)
@showbenchmark unsafe_update(y)
@showbenchmark unsafe_update2d(x2d)
@showbenchmark unsafe_update2d(y2d)
@showbenchmark unsafe_update_eachindex(x)
@showbenchmark unsafe_update_eachindex(y)

# Benchmarks of vector indexing using OffsetRanges as axes
@showbenchmark vectorlinearindexing(x, s)
@showbenchmark vectorlinearindexing(x, parent(s))
@showbenchmark vectorCartesianindexing(x2d, s, s)
@showbenchmark vectorCartesianindexing(x2d, parent(s), parent(s))
@showbenchmark nestedvectorlinearindexing(x, r, s)
@showbenchmark nestedvectorlinearindexing(x, parent(r), parent(s))