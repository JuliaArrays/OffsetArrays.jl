using Base.Test
using OffsetArrays

ambs = detect_ambiguities(Base, Core)  # in case these have ambiguities of their own
@test isempty(setdiff(detect_ambiguities(OffsetArrays, Base, Core), ambs))

# Basics
for n = 0:5
    a = OffsetArray(ones(Int,ntuple(d->1,n)), ntuple(x->x-1,n))
    @test length(linearindices(a)) == 1
    @test indices(a) == ntuple(x->x:x,n)
    @test a[1] == 1
end
a0 = reshape([3])
a = OffsetArray(a0)
@test indices(a) == ()
@test ndims(a) == 0
@test a[] == 3

y = OffsetArray(Float64, -1:1, -7:7, -128:512, -5:5, -1:1, -3:3, -2:2, -1:1)
@test indices(y) == (-1:1, -7:7, -128:512, -5:5, -1:1, -3:3, -2:2, -1:1)
y[-1,-7,-128,-5,-1,-3,-2,-1] = 14
y[-1,-7,-128,-5,-1,-3,-2,-1] += 5
@test y[-1,-7,-128,-5,-1,-3,-2,-1] == 19

r = -2:5
y = OffsetArray(r, r)
@test indices(y) == (r,)
y = OffsetArray(r, (r,))
@test indices(y) == (r,)

A0 = [1 3; 2 4]
A = OffsetArray(A0, (-1,2))                   # LinearFast
S = OffsetArray(view(A0, 1:2, 1:2), (-1,2))   # LinearSlow
@test indices(A) == indices(S) == (0:1, 3:4)
@test_throws ErrorException size(A)
@test_throws ErrorException size(A, 1)
@test A == OffsetArray(A0, 0:1, 3:4)
@test_throws DimensionMismatch OffsetArray(A0, 0:2, 3:4)
@test_throws DimensionMismatch OffsetArray(A0, 0:1, 2:4)

# Scalar indexing
@test A[0,3] == A[0,3,1] == A[1] == S[0,3] == S[0,3,1] == S[1] == 1
@test A[1,3] == A[1,3,1] == A[2] == S[1,3] == S[1,3,1] == S[2] == 2
@test A[0,4] == A[0,4,1] == A[3] == S[0,4] == S[0,4,1] == S[3] == 3
@test A[1,4] == A[1,4,1] == A[4] == S[1,4] == S[1,4,1] == S[4] == 4
@test @unsafe(A[0,3]) == @unsafe(A[0,3,1]) == @unsafe(A[1]) == @unsafe(S[0,3]) == @unsafe(S[0,3,1]) == @unsafe(S[1]) == 1
@test @unsafe(A[1,3]) == @unsafe(A[1,3,1]) == @unsafe(A[2]) == @unsafe(S[1,3]) == @unsafe(S[1,3,1]) == @unsafe(S[2]) == 2
@test @unsafe(A[0,4]) == @unsafe(A[0,4,1]) == @unsafe(A[3]) == @unsafe(S[0,4]) == @unsafe(S[0,4,1]) == @unsafe(S[3]) == 3
@test @unsafe(A[1,4]) == @unsafe(A[1,4,1]) == @unsafe(A[4]) == @unsafe(S[1,4]) == @unsafe(S[1,4,1]) == @unsafe(S[4]) == 4
@test_throws BoundsError A[1,1]
@test_throws BoundsError S[1,1]
@test_throws BoundsError A[0,3,2]
@test_throws BoundsError A[0,3,0]
Ac = copy(A)
Ac[0,3] = 10
@test Ac[0,3] == 10
Ac[0,3,1] = 11
@test Ac[0,3] == 11
@unsafe Ac[0,3,1] = 12
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
@test @unsafe(A[CartesianIndex((0,3))]) == @unsafe(S[CartesianIndex((0,3))]) == 1
@test @unsafe(A[CartesianIndex((0,3)),1]) == @unsafe(S[CartesianIndex((0,3)),1]) == 1
@test_throws BoundsError A[CartesianIndex(1,1)]
@test_throws BoundsError A[CartesianIndex(1,1),0]
@test_throws BoundsError A[CartesianIndex(1,1),2]
@test_throws BoundsError S[CartesianIndex(1,1)]
@test_throws BoundsError S[CartesianIndex(1,1),0]
@test_throws BoundsError S[CartesianIndex(1,1),2]
@test eachindex(A) == 1:4
@test eachindex(S) == CartesianRange((0:1,3:4))

# view
S = view(A, :, 3)
@test S == OffsetArray([1,2], (A.offsets[1],))
@test S[0] == 1
@test S[1] == 2
@test_throws BoundsError S[2]
@test indices(S) === (0:1,)
S = view(A, 0, :)
@test S == OffsetArray([1,3], (A.offsets[2],))
@test S[3] == 1
@test S[4] == 3
@test_throws BoundsError S[1]
@test indices(S) === (3:4,)
S = view(A, 0:0, 4)
@test S == [3]
@test S[1] == 3
@test_throws BoundsError S[0]
@test indices(S) === (Base.OneTo(1),)
S = view(A, 1, 3:4)
@test S == [2,4]
@test S[1] == 2
@test S[2] == 4
@test_throws BoundsError S[3]
@test indices(S) === (Base.OneTo(2),)
S = view(A, :, :)
@test S == A
@test S[0,3] == S[1] == 1
@test S[1,3] == S[2] == 2
@test S[0,4] == S[3] == 3
@test S[1,4] == S[4] == 4
@test_throws BoundsError S[1,1]
@test indices(S) === (0:1, 3:4)

# iteration
for (a,d) in zip(A, A0)
    @test a == d
end

# show
io = IOBuffer()
show(io, A)
str = takebuf_string(io)
@test str == "[1 3; 2 4]"
show(io, S)
str = takebuf_string(io)
@test str == "[1 3; 2 4]"
show(io, MIME("text/plain"), A)
strs = split(strip(takebuf_string(io)), '\n')
@test strs[2] == " 1  3"
@test strs[3] == " 2  4"
v = OffsetArray(rand(3), (-2,))
show(io, v)
str = takebuf_string(io)
show(io, parent(v))
@test str == takebuf_string(io)
function cmp_showf(printfunc, io, A)
    ioc = IOContext(io, limit=true, compact=true)
    printfunc(ioc, A)
    str1 = takebuf_string(io)
    printfunc(ioc, parent(A))
    str2 = takebuf_string(io)
    @test str1 == str2
end
cmp_showf(Base.print_matrix, io, OffsetArray(rand(5,5), (10,-9)))       # rows&cols fit
cmp_showf(Base.print_matrix, io, OffsetArray(rand(10^3,5), (10,-9)))    # columns fit
cmp_showf(Base.print_matrix, io, OffsetArray(rand(5,10^3), (10,-9)))    # rows fit
cmp_showf(Base.print_matrix, io, OffsetArray(rand(10^3,10^3), (10,-9))) # neither fits

# Similar
B = similar(A, Float32)
@test isa(B, OffsetArray{Float32,2})
@test indices(B) === indices(A)
B = similar(A, (3,4))
@test isa(B, Array{Int,2})
@test size(B) == (3,4)
@test indices(B) === (Base.OneTo(3), Base.OneTo(4))
B = similar(A, (-3:3,1:4))
@test isa(B, OffsetArray{Int,2})
@test indices(B) === (-3:3, 1:4)
B = similar(parent(A), (-3:3,1:4))
@test isa(B, OffsetArray{Int,2})
@test indices(B) === (-3:3, 1:4)

# Reshape
B = reshape(A0, -10:-9, 9:10)
@test isa(B, OffsetArray{Int,2})
@test parent(B) === A0
@test indices(B) == (-10:-9, 9:10)
B = reshape(A, -10:-9, 9:10)
@test isa(B, OffsetArray{Int,2})
@test parent(B) === A0
@test indices(B) == (-10:-9, 9:10)
b = reshape(A, -7:-4)
@test indices(b) == (-7:-4,)
@test isa(parent(b), Vector{Int})
@test parent(b) == A0[:]
a = OffsetArray(rand(3,3,3), -1:1, 0:2, 3:5)
@test_throws ArgumentError reshape(a, Val{2})
@test_throws ArgumentError reshape(a, Val{4})

# Indexing with OffsetArray indices
i1 = OffsetArray([2,1], (-5,))
i1 = OffsetArray([2,1], -5)
b = A0[i1, 1]
@test indices(b) === (-4:-3,)
@test b[-4] == 2
@test b[-3] == 1
b = A0[1,i1]
@test indices(b) === (-4:-3,)
@test b[-4] == 3
@test b[-3] == 1
v = view(A0, i1, 1)
@test indices(v) === (-4:-3,)
v = view(A0, 1:1, i1)
@test indices(v) === (Base.OneTo(1), -4:-3)

# logical indexing
@test A[A .> 2] == [3,4]

# copy!
a = OffsetArray{Int}((-3:-1,))
fill!(a, -1)
copy!(a, (1,2))   # non-array iterables
@test a[-3] == 1
@test a[-2] == 2
@test a[-1] == -1
fill!(a, -1)
copy!(a, -2, (1,2))
@test a[-3] == -1
@test a[-2] == 1
@test a[-1] == 2
@test_throws BoundsError copy!(a, 1, (1,2))
fill!(a, -1)
copy!(a, -2, (1,2,3), 2)
@test a[-3] == -1
@test a[-2] == 2
@test a[-1] == 3
@test_throws BoundsError copy!(a, -2, (1,2,3), 1)
fill!(a, -1)
copy!(a, -2, (1,2,3), 1, 2)
@test a[-3] == -1
@test a[-2] == 1
@test a[-1] == 2

b = 1:2    # copy between AbstractArrays
bo = OffsetArray(1:2, (-3,))
@test_throws BoundsError copy!(a, b)
fill!(a, -1)
copy!(a, bo)
@test a[-3] == -1
@test a[-2] == 1
@test a[-1] == 2
fill!(a, -1)
copy!(a, -2, bo)
@test a[-3] == -1
@test a[-2] == 1
@test a[-1] == 2
@test_throws BoundsError copy!(a, -4, bo)
@test_throws BoundsError copy!(a, -1, bo)
fill!(a, -1)
copy!(a, -3, b, 2)
@test a[-3] == 2
@test a[-2] == a[-1] == -1
@test_throws BoundsError copy!(a, -3, b, 1, 4)
am = OffsetArray{Int}((1:1, 7:9))  # for testing linear indexing
fill!(am, -1)
copy!(am, b)
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
cumsum!(C, A, 1)
@test parent(C) == cumsum(parent(A), 1)
@test parent(cumsum(A, 1)) == cumsum(parent(A), 1)
cumsum!(C, A, 2)
@test parent(C) == cumsum(parent(A), 2)
R = similar(A, (1:1, 6:9))
maximum!(R, A)
@test parent(R) == maximum(parent(A), 1)
R = similar(A, (-2:1, 1:1))
maximum!(R, A)
@test parent(R) == maximum(parent(A), 2)
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
@test isequal(cumsum_kbn(v), cv)
@test isequal(cumsum_kbn(v2), cv2)
@test isequal(sum_kbn(v), sum_kbn(parent(v)))

io = IOBuffer()
writedlm(io, A)
seek(io, 0)
@test readdlm(io, eltype(A)) == parent(A)

amin, amax = extrema(parent(A))
@test clamp.(A, (amax+amin)/2, amax) == OffsetArray(clamp.(parent(A), (amax+amin)/2, amax), indices(A))

@test unique(A, 1) == parent(A)
@test unique(A, 2) == parent(A)
v = OffsetArray(rand(8), (-2,))
@test sort(v) == OffsetArray(sort(parent(v)), v.offsets)
@test sortrows(A) == OffsetArray(sortrows(parent(A)), A.offsets)
@test sortcols(A) == OffsetArray(sortcols(parent(A)), A.offsets)
@test sort(A, 1) == OffsetArray(sort(parent(A), 1), A.offsets)
@test sort(A, 2) == OffsetArray(sort(parent(A), 2), A.offsets)

@test mapslices(v->sort(v), A, 1) == OffsetArray(mapslices(v->sort(v), parent(A), 1), A.offsets)
@test mapslices(v->sort(v), A, 2) == OffsetArray(mapslices(v->sort(v), parent(A), 2), A.offsets)

@test rotl90(A) == OffsetArray(rotl90(parent(A)), A.offsets[[2,1]])
@test rotr90(A) == OffsetArray(rotr90(parent(A)), A.offsets[[2,1]])
@test flipdim(A, 1) == OffsetArray(flipdim(parent(A), 1), A.offsets)
@test flipdim(A, 2) == OffsetArray(flipdim(parent(A), 2), A.offsets)

@test A+1 == OffsetArray(parent(A)+1, A.offsets)
@test 2*A == OffsetArray(2*parent(A), A.offsets)
@test A+A == OffsetArray(parent(A)+parent(A), A.offsets)
@test A.*A == OffsetArray(parent(A).*parent(A), A.offsets)

B = fill(5, 1:3, -1:1)
@test indices(B) == (1:3,-1:1)
@test all(B.==5)

# @unsafe
a = OffsetArray(zeros(7), -3:3)
unsafe_fill!(x) = @unsafe(for i in indices(x,1); x[i] = i; end)
function unsafe_sum(x)
    s = zero(eltype(x))
    @unsafe for i in indices(x,1)
        s += x[i]
    end
    s
end
unsafe_fill!(a)
for i = -3:3
    @test a[i] == i
end
@test unsafe_sum(a) == 0
