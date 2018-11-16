using OffsetArrays
using Test
using DelimitedFiles

@test isempty(detect_ambiguities(OffsetArrays, Base, Core))

# Basics
for n = 0:5
    for z in (OffsetArray(ones(Int,ntuple(d->1,n)), ntuple(x->x-1,n)),
              fill!(OffsetArray{Float64}(undef, ntuple(x->x:x, n)), 1),
              fill!(OffsetArray{Float64}(undef, ntuple(x->x:x, n)...), 1),
              fill!(OffsetArray{Float64,n}(undef, ntuple(x->x:x, n)), 1),
              fill!(OffsetArray{Float64,n}(undef, ntuple(x->x:x, n)...), 1))
        @test length(LinearIndices(z)) == 1
        @test axes(z) == ntuple(x->x:x, n)
        @test z[1] == 1
    end
end
a0 = reshape([3])
a = OffsetArray(a0)
@test axes(a) == ()
@test ndims(a) == 0
@test a[] == 3

# missing and nothing constructors
for (T, t) in ((Missing, missing), (Nothing, nothing))
    @test !isassigned(OffsetArray{Union{T,Vector{Int}}}(undef, -1:1, -1:1), -1, -1)
    @test OffsetArray{Union{T,Vector{Int}}}(t, -1:1, -1:1)[-1, -1] === t
    @test !isassigned(OffsetVector{Union{T,Vector{Int}}}(undef, -1:1), -1)
    @test OffsetVector{Union{T,Vector{Int}}}(t, -1:1)[-1] === t
end

y = OffsetArray{Float64}(undef, -1:1, -7:7, -128:512, -5:5, -1:1, -3:3, -2:2, -1:1)
@test axes(y) == (-1:1, -7:7, -128:512, -5:5, -1:1, -3:3, -2:2, -1:1)
y[-1,-7,-128,-5,-1,-3,-2,-1] = 14
y[-1,-7,-128,-5,-1,-3,-2,-1] += 5
@test y[-1,-7,-128,-5,-1,-3,-2,-1] == 19

r = -2:5
y = OffsetArray(r, r)
@test axes(y) == (r,)
y = OffsetArray(r, (r,))
@test axes(y) == (r,)

y = OffsetArray{Float32}(undef, (Base.Slice(-1:1),))
@test axes(y) === (Base.Slice(-1:1),)

A0 = [1 3; 2 4]
A = OffsetArray(A0, (-1,2))                   # IndexLinear
S = OffsetArray(view(A0, 1:2, 1:2), (-1,2))   # IndexCartesian
@test axes(A) == axes(S) == (0:1, 3:4)
@test size(A) == size(A0)
@test size(A, 1) == size(A0, 1)
@test length(A) == length(A0)
@test A == OffsetArray(A0, 0:1, 3:4)
@test_throws DimensionMismatch OffsetArray(A0, 0:2, 3:4)
@test_throws DimensionMismatch OffsetArray(A0, 0:1, 2:4)

# Scalar indexing
@test @inferred(A[0,3]) == @inferred(A[0,3,1]) == @inferred(A[1]) == @inferred(S[0,3]) == @inferred(S[0,3,1]) == @inferred(S[1]) == 1
@test A[1,3] == A[1,3,1] == A[2] == S[1,3] == S[1,3,1] == S[2] == 2
@test A[0,4] == A[0,4,1] == A[3] == S[0,4] == S[0,4,1] == S[3] == 3
@test A[1,4] == A[1,4,1] == A[4] == S[1,4] == S[1,4,1] == S[4] == 4
@test @inbounds(A[0,3]) == @inbounds(A[0,3,1]) == @inbounds(A[1]) == @inbounds(S[0,3]) == @inbounds(S[0,3,1]) == @inbounds(S[1]) == 1
@test @inbounds(A[1,3]) == @inbounds(A[1,3,1]) == @inbounds(A[2]) == @inbounds(S[1,3]) == @inbounds(S[1,3,1]) == @inbounds(S[2]) == 2
@test @inbounds(A[0,4]) == @inbounds(A[0,4,1]) == @inbounds(A[3]) == @inbounds(S[0,4]) == @inbounds(S[0,4,1]) == @inbounds(S[3]) == 3
@test @inbounds(A[1,4]) == @inbounds(A[1,4,1]) == @inbounds(A[4]) == @inbounds(S[1,4]) == @inbounds(S[1,4,1]) == @inbounds(S[4]) == 4
@test_throws BoundsError A[1,1]
@test_throws BoundsError S[1,1]
@test_throws BoundsError A[0,3,2]
@test_throws BoundsError A[0,3,0]
Ac = copy(A)
Ac[0,3] = 10
@test Ac[0,3] == 10
Ac[0,3,1] = 11
@test Ac[0,3] == 11
@inbounds Ac[0,3,1] = 12
@test Ac[0,3] == 12

# Vector indexing
@test A[:, 3] == S[:, 3] == OffsetArray([1,2], (A.offsets[1],))
@test A[:, 4] == S[:, 4] == OffsetArray([3,4], (A.offsets[1],))
@test_throws BoundsError A[:, 1]
@test_throws BoundsError S[:, 1]
@test A[0, :] == S[0, :] == OffsetArray([1,3], (A.offsets[2],))
@test A[1, :] == S[1, :] == OffsetArray([2,4], (A.offsets[2],))
@test_throws BoundsError A[2, :]
@test_throws BoundsError S[2, :]
@test A[0:1, 3] == S[0:1, 3] == [1,2]
@test A[[1,0], 3] == S[[1,0], 3] == [2,1]
@test A[0, 3:4] == S[0, 3:4] == [1,3]
@test A[1, [4,3]] == S[1, [4,3]] == [4,2]
@test A[:, :] == S[:, :] == A

# CartesianIndexing
@test A[CartesianIndex((0,3))] == S[CartesianIndex((0,3))] == 1
@test A[CartesianIndex((0,3)),1] == S[CartesianIndex((0,3)),1] == 1
@test @inbounds(A[CartesianIndex((0,3))]) == @inbounds(S[CartesianIndex((0,3))]) == 1
@test @inbounds(A[CartesianIndex((0,3)),1]) == @inbounds(S[CartesianIndex((0,3)),1]) == 1
@test_throws BoundsError A[CartesianIndex(1,1)]
@test_throws BoundsError A[CartesianIndex(1,1),0]
@test_throws BoundsError A[CartesianIndex(1,1),2]
@test_throws BoundsError S[CartesianIndex(1,1)]
@test_throws BoundsError S[CartesianIndex(1,1),0]
@test_throws BoundsError S[CartesianIndex(1,1),2]
@test eachindex(A) == 1:4
@test eachindex(S) == CartesianIndices(Base.Slice.((0:1,3:4)))

# view
S = view(A, :, 3)
@test S == OffsetArray([1,2], (A.offsets[1],))
@test S[0] == 1
@test S[1] == 2
@test_throws BoundsError S[2]
@test axes(S) === (Base.Slice(0:1),)
S = view(A, 0, :)
@test S == OffsetArray([1,3], (A.offsets[2],))
@test S[3] == 1
@test S[4] == 3
@test_throws BoundsError S[1]
@test axes(S) === (Base.Slice(3:4),)
S = view(A, 0:0, 4)
@test S == [3]
@test S[1] == 3
@test_throws BoundsError S[0]
@test axes(S) === (Base.OneTo(1),)
S = view(A, 1, 3:4)
@test S == [2,4]
@test S[1] == 2
@test S[2] == 4
@test_throws BoundsError S[3]
@test axes(S) === (Base.OneTo(2),)
S = view(A, :, :)
@test S == A
@test S[0,3] == S[1] == 1
@test S[1,3] == S[2] == 2
@test S[0,4] == S[3] == 3
@test S[1,4] == S[4] == 4
@test_throws BoundsError S[1,1]
@test axes(S) === Base.Slice.((0:1, 3:4))

# iteration
let a
    for (a,d) in zip(A, A0)
        @test a == d
    end
end

# show
@test sprint(show, A) == "[1 3; 2 4]"
@test sprint(show, S) == "[1 3; 2 4]"
strs = split(strip(sprint(show, MIME("text/plain"), A)), '\n')
@test strs[2] == " 1  3"
@test strs[3] == " 2  4"
v = OffsetArray(rand(3), (-2,))
@test sprint(show, v) == sprint(show, parent(v))
io = IOBuffer()
function cmp_showf(printfunc, io, A)
    ioc = IOContext(io, :limit=>true, :compact=>true)
    printfunc(ioc, A)
    str1 = String(take!(io))
    printfunc(ioc, parent(A))
    str2 = String(take!(io))
    @test str1 == str2
end
cmp_showf(Base.print_matrix, io, OffsetArray(rand(5,5), (10,-9)))       # rows&cols fit
cmp_showf(Base.print_matrix, io, OffsetArray(rand(10^3,5), (10,-9)))    # columns fit
cmp_showf(Base.print_matrix, io, OffsetArray(rand(5,10^3), (10,-9)))    # rows fit
cmp_showf(Base.print_matrix, io, OffsetArray(rand(10^3,10^3), (10,-9))) # neither fits

# Similar
B = similar(A, Float32)
@test isa(B, OffsetArray{Float32,2})
@test axes(B) === axes(A)
B = similar(A, (3,4))
@test isa(B, Array{Int,2})
@test size(B) == (3,4)
@test axes(B) === (Base.OneTo(3), Base.OneTo(4))
B = similar(A, (-3:3,1:4))
@test isa(B, OffsetArray{Int,2})
@test axes(B) === Base.Slice.((-3:3, 1:4))
B = similar(parent(A), (-3:3,1:4))
@test isa(B, OffsetArray{Int,2})
@test axes(B) === Base.Slice.((-3:3, 1:4))
@test isa([x for x in [1,2,3]], Vector{Int})
@test similar(Array{Int}, (0:0, 0:0)) isa OffsetArray{Int, 2}
@test similar(Array{Int}, (1, 1)) isa Matrix{Int}
@test similar(Array{Int}, (Base.OneTo(1), Base.OneTo(1))) isa Matrix{Int}

# Reshape
B = reshape(A0, -10:-9, 9:10)
@test isa(B, OffsetArray{Int,2})
@test parent(B) === A0
@test axes(B) == Base.Slice.((-10:-9, 9:10))
B = reshape(A, -10:-9, 9:10)
@test isa(B, OffsetArray{Int,2})
@test pointer(parent(B)) === pointer(A0)
@test axes(B) == Base.Slice.((-10:-9, 9:10))
b = reshape(A, -7:-4)
@test axes(b) == (Base.Slice(-7:-4),)
@test isa(parent(b), Vector{Int})
@test pointer(parent(b)) === pointer(parent(A))
@test parent(b) == A0[:]
a = OffsetArray(rand(3,3,3), -1:1, 0:2, 3:5)
# Offset axes are required for reshape(::OffsetArray, ::Val) support
b = reshape(a, Val(2))
@test isa(b, OffsetArray{Float64,2})
@test pointer(parent(b)) === pointer(parent(a))
@test axes(b) == Base.Slice.((-1:1, 1:9))
b = reshape(a, Val(4))
@test isa(b, OffsetArray{Float64,4})
@test pointer(parent(b)) === pointer(parent(a))
@test axes(b) == (axes(a)..., Base.Slice(1:1))

# Indexing with OffsetArray axes
i1 = OffsetArray([2,1], (-5,))
i1 = OffsetArray([2,1], -5)
b = A0[i1, 1]
@test axes(b) === (Base.Slice(-4:-3),)
@test b[-4] == 2
@test b[-3] == 1
b = A0[1,i1]
@test axes(b) === (Base.Slice(-4:-3),)
@test b[-4] == 3
@test b[-3] == 1
v = view(A0, i1, 1)
@test axes(v) === (Base.Slice(-4:-3),)
v = view(A0, 1:1, i1)
@test axes(v) === (Base.OneTo(1), Base.Slice(-4:-3))

# logical indexing
@test A[A .> 2] == [3,4]

# copyto!
a = OffsetArray{Int}(undef, (-3:-1,))
fill!(a, -1)
copyto!(a, (1,2))   # non-array iterables
@test a[-3] == 1
@test a[-2] == 2
@test a[-1] == -1
fill!(a, -1)
copyto!(a, -2, (1,2))
@test a[-3] == -1
@test a[-2] == 1
@test a[-1] == 2
@test_throws BoundsError copyto!(a, 1, (1,2))
fill!(a, -1)
copyto!(a, -2, (1,2,3), 2)
@test a[-3] == -1
@test a[-2] == 2
@test a[-1] == 3
@test_throws BoundsError copyto!(a, -2, (1,2,3), 1)
fill!(a, -1)
copyto!(a, -2, (1,2,3), 1, 2)
@test a[-3] == -1
@test a[-2] == 1
@test a[-1] == 2

b = 1:2    # copy between AbstractArrays
bo = OffsetArray(1:2, (-3,))
@test_throws BoundsError copyto!(a, b)
fill!(a, -1)
copyto!(a, bo)
@test a[-3] == -1
@test a[-2] == 1
@test a[-1] == 2
fill!(a, -1)
copyto!(a, -2, bo)
@test a[-3] == -1
@test a[-2] == 1
@test a[-1] == 2
@test_throws BoundsError copyto!(a, -4, bo)
@test_throws BoundsError copyto!(a, -1, bo)
fill!(a, -1)
copyto!(a, -3, b, 2)
@test a[-3] == 2
@test a[-2] == a[-1] == -1
@test_throws BoundsError copyto!(a, -3, b, 1, 4)
am = OffsetArray{Int}(undef, (1:1, 7:9))  # for testing linear indexing
fill!(am, -1)
copyto!(am, b)
@test am[1] == 1
@test am[2] == 2
@test am[3] == -1
@test am[1,7] == 1
@test am[1,8] == 2
@test am[1,9] == -1

dest = similar(am)
map!(+, dest, am, am)
@test dest[1,7] == 2
@test dest[1,8] == 4
@test dest[1,9] == -2

A = OffsetArray(rand(4,4), (-3,5))
@test maximum(A) == maximum(parent(A))
@test minimum(A) == minimum(parent(A))
@test extrema(A) == extrema(parent(A))
C = similar(A)
cumsum!(C, A, dims = 1)
@test parent(C) == cumsum(parent(A), dims = 1)
@test parent(cumsum(A, dims = 1)) == cumsum(parent(A), dims = 1)
cumsum!(C, A, dims = 2)
@test parent(C) == cumsum(parent(A), dims = 2)
R = similar(A, (1:1, 6:9))
maximum!(R, A)
@test parent(R) == maximum(parent(A), dims = 1)
R = similar(A, (-2:1, 1:1))
maximum!(R, A)
@test parent(R) == maximum(parent(A), dims = 2)
amin, iamin = findmin(A)
pmin, ipmin = findmin(parent(A))
@test amin == pmin
@test A[iamin] == amin
@test amin == parent(A)[ipmin]
amax, iamax = findmax(A)
pmax, ipmax = findmax(parent(A))
@test amax == pmax
@test A[iamax] == amax
@test amax == parent(A)[ipmax]

v  = OffsetArray([1,1e100,1,-1e100], (-3,))*1000
v2 = OffsetArray([1,-1e100,1,1e100], (5,))*1000
@test isa(v, OffsetArray)
cv  = OffsetArray([1,1e100,1e100,2], (-3,))*1000
cv2 = OffsetArray([1,-1e100,-1e100,2], (5,))*1000
# @test isequal(cumsum_kbn(v), cv)
# @test isequal(cumsum_kbn(v2), cv2)
# @test isequal(sum_kbn(v), sum_kbn(parent(v)))

io = IOBuffer()
writedlm(io, A)
seek(io, 0)
@test readdlm(io, eltype(A)) == parent(A)

amin, amax = extrema(parent(A))
@test clamp.(A, (amax+amin)/2, amax) == OffsetArray(clamp.(parent(A), (amax+amin)/2, amax), axes(A))

@test unique(A, dims=1) == OffsetArray(parent(A), 0, first(axes(A, 2)) - 1)
@test unique(A, dims=2) == OffsetArray(parent(A), first(axes(A, 1)) - 1, 0)
v = OffsetArray(rand(8), (-2,))
@test sort(v) == OffsetArray(sort(parent(v)), v.offsets)
@test sortslices(A; dims=1) == OffsetArray(sortslices(parent(A); dims=1), A.offsets)
@test sortslices(A; dims=2) == OffsetArray(sortslices(parent(A); dims=2), A.offsets)
@test sort(A, dims = 1) == OffsetArray(sort(parent(A), dims = 1), A.offsets)
@test sort(A, dims = 2) == OffsetArray(sort(parent(A), dims = 2), A.offsets)

@test mapslices(v->sort(v), A, dims = 1) == OffsetArray(mapslices(v->sort(v), parent(A), dims = 1), A.offsets)
@test mapslices(v->sort(v), A, dims = 2) == OffsetArray(mapslices(v->sort(v), parent(A), dims = 2), A.offsets)

@test rotl90(A) == OffsetArray(rotl90(parent(A)), A.offsets[[2,1]])
@test rotr90(A) == OffsetArray(rotr90(parent(A)), A.offsets[[2,1]])
@test reverse(A, dims = 1) == OffsetArray(reverse(parent(A), dims = 1), A.offsets)
@test reverse(A, dims = 2) == OffsetArray(reverse(parent(A), dims = 2), A.offsets)

@test A.+1 == OffsetArray(parent(A).+1, A.offsets)
@test 2*A == OffsetArray(2*parent(A), A.offsets)
@test A+A == OffsetArray(parent(A)+parent(A), A.offsets)
@test A.*A == OffsetArray(parent(A).*parent(A), A.offsets)

B = fill(5, 1:3, -1:1)
@test axes(B) == (1:3,-1:1)
@test all(B.==5)

# @inbounds
a = OffsetArray(zeros(7), -3:3)
unsafe_fill!(x) = @inbounds(for i in axes(x,1); x[i] = i; end)
function unsafe_sum(x)
    s = zero(eltype(x))
    @inbounds for i in axes(x,1)
        s += x[i]
    end
    s
end
unsafe_fill!(a)
for i = -3:3
    @test a[i] == i
end
@test unsafe_sum(a) == 0

a = OffsetArray([1 2; 3 4], -1:0, 5:6)
@test summary(a) == "OffsetArray(::Array{$(Int),2}, -1:0, 5:6) with eltype $(Int) with indices -1:0×5:6"
@test summary(view(a, :, 5)) == "view(OffsetArray(::Array{$(Int),2}, -1:0, 5:6), :, 5) with eltype $(Int) with indices -1:0"
a = OffsetArray(reshape([1]))
@test summary(a) == "0-dimensional OffsetArray(::Array{$(Int),0}) with eltype $(Int)"

@testset "OffsetVector constructors" begin
    local v = rand(5)
    @test OffsetVector(v, -2) == OffsetArray(v, -2)
    @test OffsetVector(v, -2:2) == OffsetArray(v, -2:2)
    @test typeof(OffsetVector{Float64}(undef, -2:2)) == typeof(OffsetArray{Float64}(undef, -2:2))
end
