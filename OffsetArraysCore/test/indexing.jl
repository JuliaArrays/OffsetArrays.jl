module IndexingTests

using OffsetArraysCore
using OffsetArraysCore: OffsetArray, IdOffsetRange, no_offset_view, IdentityUnitRange
using StaticArrays
using FillArrays
using Test

no_offset_axes(x, d) = no_offset_view(axes(x, d))
no_offset_axes(x) = map(no_offset_view, axes(x))

function same_value(r1, r2)
    length(r1) == length(r2) || return false
    for (v1, v2) in zip(r1, r2)
        v1 == v2 || return false
    end
    return true
end

include("customranges.jl")

@testset "Traits" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))                   # IndexLinear
    S = OffsetArray(view(A0, 1:2, 1:2), (-1,2))   # IndexCartesian
    @test axes(A) === axes(S)
    @test no_offset_axes(A) == no_offset_axes(S) == (0:1, 3:4)
    @test axes(A, 1) === IdOffsetRange(Base.OneTo(2), -1)
    @test size(A) == size(A0)
    @test size(A, 1) == size(A0, 1)
    @test length(A) == length(A0)
    @test A == OffsetArray(A0, 0:1, 3:4)
    @test_throws DimensionMismatch OffsetArray(A0, 0:2, 3:4)
    @test_throws DimensionMismatch OffsetArray(A0, 0:1, 2:4)
    @test eachindex(IndexLinear(), A) == eachindex(IndexLinear(), parent(A))
    @test eachindex(IndexCartesian(), A) == CartesianIndices(A) == CartesianIndices(axes(A))
    @test eachindex(S) == eachindex(IndexCartesian(), S) == CartesianIndices(S)
    @test eachindex(IndexLinear(), S) == eachindex(IndexLinear(), A0)
    A = OffsetArray(ones(2), 5:6)
    @test eachindex(IndexLinear(), A) === axes(A, 1)

    A = OffsetArray(big(1):big(2), 1)
    B = OffsetArray(1:2, 1)
    @test CartesianIndices(A) == CartesianIndices(B)
    @test LinearIndices(A) == LinearIndices(B)
    @test eachindex(A) == eachindex(B)
end

@testset "Scalar indexing" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))
    S = OffsetArray(view(A0, 1:2, 1:2), (-1,2))

    @test @inferred(A[0,3]) == @inferred(A[0,3,1]) == @inferred(A[1]) == @inferred(S[0,3]) == @inferred(S[0,3,1]) == @inferred(S[1]) == 1
    @test A[1,3] == A[1,3,1] == A[2] == S[1,3] == S[1,3,1] == S[2] == 2
    @test A[0,4] == A[0,4,1] == A[3] == S[0,4] == S[0,4,1] == S[3] == 3
    @test A[1,4] == A[1,4,1] == A[4] == S[1,4] == S[1,4,1] == S[4] == 4
    @test @inbounds(A[0,3]) == @inbounds(A[0,3,1]) == @inbounds(A[1]) == @inbounds(S[0,3]) == @inbounds(S[0,3,1]) == @inbounds(S[1]) == 1
    @test @inbounds(A[1,3]) == @inbounds(A[1,3,1]) == @inbounds(A[2]) == @inbounds(S[1,3]) == @inbounds(S[1,3,1]) == @inbounds(S[2]) == 2
    @test @inbounds(A[0,4]) == @inbounds(A[0,4,1]) == @inbounds(A[3]) == @inbounds(S[0,4]) == @inbounds(S[0,4,1]) == @inbounds(S[3]) == 3
    @test @inbounds(A[1,4]) == @inbounds(A[1,4,1]) == @inbounds(A[4]) == @inbounds(S[1,4]) == @inbounds(S[1,4,1]) == @inbounds(S[4]) == 4
    @test_throws BoundsError(A, (1,1)) A[1,1]
    @test_throws BoundsError(A, (1,1)) A[1,1] = 4
    @test_throws BoundsError(S, (1,1)) S[1,1]
    @test_throws BoundsError(S, (1,1)) S[1,1] = 4
    @test_throws BoundsError(A, (0,3,2)) A[0,3,2]
    @test_throws BoundsError(A, (0,3,2)) A[0,3,2] = 4
    @test_throws BoundsError(A, (0,3,0)) A[0,3,0]
    @test_throws BoundsError(A, (0,3,0)) A[0,3,0] = 4
    Ac = copy(A)
    Ac[0,3] = 10
    @test Ac[0,3] == 10
    Ac[0,3,1] = 11
    @test Ac[0,3] == 11
    @inbounds Ac[0,3,1] = 12
    @test Ac[0,3] == 12

    y = OffsetArray{Float64}(undef, -1:1, -7:7, -3:-1, -5:5, -1:1, -3:3, -2:2, -1:1)
    y[-1,-7,-3,-5,-1,-3,-2,-1] = 14
    y[-1,-7,-3,-5,-1,-3,-2,-1] += 5
    @test y[-1,-7,-3,-5,-1,-3,-2,-1] == 19

    @testset "setindex!" begin
        A = OffsetArray(ones(2,2), 1:2, 1:2)
        @test setindex!(A, 2, 1, 1) === A
        @test A[1,1] == 2
        @test setindex!(A, 2, 1) === A
        @test A[1] == 2

        v = OffsetArray(ones(3), 4:6)
        @test setindex!(A, 2, 4) === A
        @test A[4] == 2
    end

    @testset "Zero-index indexing (#194)" begin
        @test OffsetArray([6], 2:2)[] == 6
        @test OffsetArray(fill(6, 1, 1), 2:2, 3:3)[] == 6
        @test OffsetArray(fill(6))[] == 6
        @test_throws BoundsError OffsetArray([6,7], 2:3)[]
        @test_throws BoundsError OffsetArray([6 7], 2:2, 2:3)[]
        @test_throws BoundsError OffsetArray([], 2:1)[]
    end
end


_comp(x::Integer, y::Integer) = x == y
_comp(x::Any, y::Any) = isapprox(Real(x), Real(y), atol = 1e-14, rtol = 1e-8)
function test_indexing_axes_and_vals(r1, r2)
    r12 = r1[r2]
    if axes(r12, 1) != axes(r2, 1)
        @show r1 r2 r12 axes(r12, 1) axes(r2, 1)
    end
    @test axes(r12, 1) == axes(r2, 1)

    if axes(r12, 1) == axes(r2, 1)
        res1 = try
            _comp(first(r12), r1[first(r2)])
        catch
            @show r1 r2
            rethrow()
        end
        res2 = try
            _comp(last(r12), r1[last(r2)])
        catch
            @show r1 r2
            rethrow()
        end
        if !(res1 & res2)
            @show r1 r2
        end
        @test res1
        @test res2
        for i in eachindex(r2)
            @test _comp(r12[i], r1[r2[i]])
        end
    end
end

@testset "Vector indexing" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))
    S = OffsetArray(view(A0, 1:2, 1:2), (-1,2))

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

    # Indexing a nD OffsetArray with n colons preserves the type
    r1 = OffsetArray(IdentityUnitRange(100:1000), 3)
    @test r1[:] === r1

    # In general with more colons than dimensions,
    # the type might not be preserved but the values and the leading axes should be
    r2 = r1[:,:]
    @test axes(r2, 1) == axes(r1, 1)
    @test same_value(r1, r2)

    s = @SVector[i for i in 1:3]
    so = OffsetArray(s, 3)
    @test so[:] === so

    a = Ones(3, 2, 1)
    ao = OffsetArray(a, axes(a))
    @test ao[:,:,:] === ao
    @test same_value(ao[:], ao)
    @test same_value(ao[:,:], ao)

    # Indexing an nD OffsetArray with one Colon preserves only the values.
    # This uses linear indexing
    a = OffsetArray(ones(2,2), (2:3, 2:3))
    b = a[:]
    @test same_value(a, b)

    vals = (1,2,3,4,5,6,7,8)
    s = SArray{Tuple{2,2,2},Int,3,8}(vals)
    so = OffsetArray(s, axes(s));
    so2 = so[:]
    @test same_value(so2, s)

    # Test r1[inds] for various combinations of types

    # AbstractArrays with 1-based indices
    indslist1 = Any[
            OffsetArray(5:8, 0),
            # This currently errors for IdentityUnitRange
            # see https://github.com/JuliaLang/julia/issues/39997
            # OffsetArray(big(5):big(80), 0),
            OffsetArray(5:2:9, 0),
            OffsetArray(9:-2:5, 0),
            OffsetArray(IdentityUnitRange(5:8), -4),
            OffsetArray(IdOffsetRange(5:8), 0),
            ]

    # AbstractRanges with 1-based indices
    indslist2 = Any[
            5:8,
            # This currently errors for IdentityUnitRange
            # see https://github.com/JuliaLang/julia/issues/39997
            # big(5):big(80),
            5:2:9,
            9:-2:5,
            IdOffsetRange(5:8),
            IdOffsetRange(ZeroBasedUnitRange(4:7), 1),
            ]

    for r1 in Any[
        # AbstractArrays
        collect(1:100),
        OffsetArray(collect(-1:100), -1),
        OffsetArray(collect(reshape(1:400, 20, 20)), -10, -10),
        OffsetArray(collect(reshape(1:21^2, 21, 21)), -10, -10),

        # OffsetRanges
        OffsetArray(10:1000, 0), # 1-based index
        OffsetArray(UnitRange(10.0, 1000.0), 0), # 1-based index
        OffsetArray(10:3:1000, 3), # offset index
        OffsetArray(10.0:3:1000.0, 0), # 1-based index
        OffsetArray(10.0:3:1000.0, 3), # offset index
        OffsetArray(IdOffsetRange(10:1000, 1), -1), # 1-based index
        OffsetArray(IdOffsetRange(10:1000, 1), 3), # offset index
        OffsetArray(IdOffsetRange(IdOffsetRange(10:1000, -4), 1), 3), # 1-based index
        OffsetArray(IdOffsetRange(IdOffsetRange(10:1000, -1), 1), 3), # offset index

        # AbstractRanges
        Base.OneTo(1000),
        CustomRange(Base.OneTo(1000)),
        Base.Slice(Base.OneTo(1000)),
        1:1000,
        UnitRange(1.0, 1000.0),
        1:3:1000,
        1000:-3:1,
        1.0:3.0:1000.0,
        StepRangeLen(Float64(1), Float64(1000), 1000),
        LinRange(1, 1000, 1000),
        Base.Slice(Base.OneTo(1000)), # 1-based index
        IdentityUnitRange(Base.OneTo(1000)), # 1-based index
        IdOffsetRange(Base.OneTo(1000)), # 1-based index
        IdentityUnitRange(2:1000), # offset index
        IdOffsetRange(ZeroBasedUnitRange(1:1000), 1), # 1-based index
        IdOffsetRange(ZeroBasedUnitRange(1:1000), 2), # offset index
        ZeroBasedUnitRange(1:1000), # offset range
        ZeroBasedRange(1:1000), # offset range
        ZeroBasedRange(1:1:1000), # offset range
        CustomRange(ZeroBasedRange(1:1:1000)), # offset range
        ]

        # AbstractArrays with 1-based indices
        for r2 in indslist1
            test_indexing_axes_and_vals(r1, r2)
            test_indexing_axes_and_vals(r1, collect(r2))
        end

        # AbstractRanges with 1-based indices
        for r2 in indslist2
            if !(!(r1 isa AbstractUnitRange) && r1 isa AbstractRange && r2 isa OffsetVector)
                test_indexing_axes_and_vals(r1, r2)
                test_indexing_axes_and_vals(r1, collect(r2))
            end

            if r1 isa AbstractRange && !(r1 isa CustomRange) && axes(r2, 1) isa Base.OneTo
                @test r1[r2] isa AbstractRange
            end
        end
    end

    # Indexing with IdentityUnitRange(::Base.OneTo) or Base.Slice(::OneTo) is special.
    # This is because axes(::IdentityUnitRange{<:Base.OneTo}, 1) isa Base.OneTo, and not an IdentityUnitRange.
    # These therefore may pass through no_offset_view unchanged.
    # This had led to a stack-overflow in indexing, as getindex was using no_offset_view.
    # Issue 209
    for r1 in Any[
        # This set of tests is for ranges r1 that have 1-based indices
        UnitRange(1.0, 99.0),
        1:99,
        Base.OneTo(99),
        1:1:99,
        99:-1:1,
        1.0:1.0:99.0,
        StepRangeLen(Float64(1), Float64(99), 99),
        LinRange(1, 99, 99),
        Base.Slice(Base.OneTo(99)),
        IdentityUnitRange(Base.OneTo(99)),
        IdOffsetRange(Base.OneTo(99)),
        ]

        for r2 in Any[
            IdentityUnitRange(Base.OneTo(3)),
            Base.Slice(Base.OneTo(3)),
            IdOffsetRange(Base.OneTo(3)),
            ]

            test_indexing_axes_and_vals(r1, r2)
            test_indexing_axes_and_vals(r1, collect(r2))
            if axes(r2, 1) isa Base.OneTo
                @test r1[r2] isa AbstractRange
            end
        end
    end
end

@testset "Vector indexing with offset ranges" begin
    r = OffsetArray(8:10, -1:1)
    r1 = r[0:1]
    @test r1 === 9:10
    r1 = (8:10)[OffsetArray(1:2, -5:-4)]
    @test no_offset_axes(r1) == (-5:-4,)
    @test no_offset_view(r1) == 8:9
    r1 = OffsetArray(8:10, -1:1)[OffsetArray(0:1, -5:-4)]
    @test no_offset_axes(r1) == (-5:-4,)
    @test no_offset_view(r1) == 9:10

    a = OffsetVector(3:4, 10:11)
    ax = IdOffsetRange(5:6, 5)
    @test axes(a[ax]) == axes(ax)
    for i in axes(ax,1)
        @test a[ax[i]] == a[ax][i]
    end

    ax = IdentityUnitRange(10:11)
    @test axes(a[ax]) == axes(ax)
    for i in axes(ax,1)
        @test a[ax[i]] == a[ax][i]
    end

    # AbstractArrays with offset axes
    indslist1 = Any[OffsetArray(5:9, 40),
            OffsetArray(5:2:9, 40),
            OffsetArray(9:-2:5, 40),
            OffsetArray(IdentityUnitRange(5:8), 2),
            OffsetArray(IdOffsetRange(5:8, 1), 3),
            OffsetArray(IdOffsetRange(IdOffsetRange(5:8, 4), 1), 3),
            OffsetArray(IdOffsetRange(IdentityUnitRange(5:8), 1), 3),
            OffsetArray(IdentityUnitRange(IdOffsetRange(5:8, 1)), 3),
            ]

    # AbstractRanges with offset axes
    indslist2 = Any[IdOffsetRange(5:8, 1),
            IdOffsetRange(Base.OneTo(3), 4),
            IdOffsetRange(IdOffsetRange(5:8, 2), 1),
            IdOffsetRange(IdOffsetRange(IdOffsetRange(5:8, -1), 2), 1),
            IdOffsetRange(IdentityUnitRange(15:20), -2),
            ]

    indslist3 = Any[IdentityUnitRange(5:8),
            IdentityUnitRange(IdOffsetRange(1:4, 5)),
            ZeroBasedUnitRange(5:8),
            ZeroBasedRange(5:8),
            ZeroBasedRange(5:2:9),
            ZeroBasedRange(9:-2:5),
            ]

    for r1 in Any[
        # AbstractArrays
        collect(1:100),
        OffsetArray(collect(-1:100), -1),
        OffsetArray(collect(reshape(1:400, 20, 20)), -10, -10),
        OffsetArray(collect(reshape(1:21^2, 21, 21)), -10, -10),

        # OffsetRanges
        OffsetArray(10:1000, 0), # 1-based index
        OffsetArray(UnitRange(10.0, 1000.0), 0), # 1-based index
        OffsetArray(10:1000, 3), # offset index
        OffsetArray(10:3:1000, 0), # 1-based index
        OffsetArray(10:3:1000, 3), # offset index
        OffsetArray(10.0:3:1000.0, 0), # 1-based index
        OffsetArray(10.0:3:1000.0, 3), # offset index
        OffsetArray(IdOffsetRange(10:1000, -3), 3), # 1-based index
        OffsetArray(IdOffsetRange(10:1000, 1), 3), # offset index
        OffsetArray(IdOffsetRange(IdOffsetRange(10:1000, -4), 1), 3), # 1-based index
        OffsetArray(IdOffsetRange(IdOffsetRange(10:1000, -1), 1), 3), # offset index

        # AbstractRanges
        Base.OneTo(1000),
        Base.Slice(Base.OneTo(1000)),
        CustomRange(Base.OneTo(1000)),
        1:1000,
        UnitRange(1.0, 1000.0),
        1:2:2000,
        2000:-1:1,
        1.0:2.0:2000.0,
        StepRangeLen(Float64(1), Float64(1000), 1000),
        LinRange(1.0, 2000.0, 2000),
        Base.Slice(Base.OneTo(1000)), # 1-based index
        IdOffsetRange(Base.OneTo(1000)), # 1-based index
        IdOffsetRange(1:1000, 0), # 1-based index
        IdOffsetRange(Base.OneTo(1000), 4), # offset index
        IdOffsetRange(1:1000, 4), # offset index
        IdOffsetRange(ZeroBasedUnitRange(1:1000), 1), # 1-based index
        IdOffsetRange(ZeroBasedUnitRange(1:1000), 2), # offset index
        IdentityUnitRange(ZeroBasedUnitRange(1:1000)), # 1-based index
        IdentityUnitRange(5:1000), # offset index
        ZeroBasedUnitRange(1:1000), # offset index
        ZeroBasedRange(1:1000), # offset index
        ZeroBasedRange(1:1:1000), # offset index
        ZeroBasedUnitRange(IdentityUnitRange(1:1000)), # offset index
        CustomRange(ZeroBasedUnitRange(IdentityUnitRange(1:1000))), # offset index
        ]

        # AbstractArrays with offset axes
        for r2 in indslist1
            test_indexing_axes_and_vals(r1, r2)
            r2_dense = OffsetArray(collect(r2), axes(r2))
            test_indexing_axes_and_vals(r1, r2_dense)
        end

        # AbstractRanges with offset axes
        for r2 in indslist2

            test_indexing_axes_and_vals(r1, r2)
            r2_dense = OffsetArray(collect(r2), axes(r2))
            test_indexing_axes_and_vals(r1, r2_dense)

            # This might not hold for all ranges, but holds for the known ones being tested here
            if r1 isa AbstractUnitRange{<:Integer} && r2 isa AbstractUnitRange{<:Integer}
                @test r1[r2] isa AbstractUnitRange{<:Integer}
            end
        end
    end
end

@testset "LinearIndexing" begin
    r = OffsetArray(ZeroBasedRange(3:4), 1)
    @test LinearIndices(r) == axes(r,1)
    r = OffsetArray(ZeroBasedRange(3:4), 2)
    @test LinearIndices(r) == axes(r,1)
end

@testset "CartesianIndexing" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))
    S = OffsetArray(view(A0, 1:2, 1:2), (-1,2))

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
    @test eachindex(S) == CartesianIndices(IdentityUnitRange.((0:1,3:4)))
end

@testset "IdentityUnitRange indexing" begin
    # 155
    a = OffsetVector(3:4, 2:3)
    ax = IdentityUnitRange(2:3)
    @test a[ax[2]] == a[ax][2]

    s = -2:2:4
    r = 5:8
    y = OffsetArray(s, r)
    @test no_offset_axes(y) == (r,)
    @test step(y) == step(s)

    a = OffsetVector(3:4, 10:11)
    ax = IdOffsetRange(5:6, 5)
    @test axes(a[ax]) == axes(ax)
    for i in axes(ax,1)
        @test a[ax[i]] == a[ax][i]
    end

    ax = IdentityUnitRange(10:11)
    @test axes(a[ax]) == axes(ax)
    for i in axes(ax,1)
        @test a[ax[i]] == a[ax][i]
    end
end

@testset "Indexing with OffsetArray axes" begin
    A0 = [1 3; 2 4]

    i1 = OffsetArray([2,1], (-5,))
    i1 = OffsetArray([2,1], -5)
    b = A0[i1, 1]
    @test axes(b) == (IdentityUnitRange(-4:-3),)
    @test b[-4] == 2
    @test b[-3] == 1
    b = A0[1,i1]
    @test axes(b) == (IdentityUnitRange(-4:-3),)
    @test b[-4] == 3
    @test b[-3] == 1
    v = view(A0, i1, 1)
    @test axes(v) == (IdentityUnitRange(-4:-3),)
    v = view(A0, 1:1, i1)
    @test axes(v) == (Base.OneTo(1), IdentityUnitRange(-4:-3))

    for r in (1:10, 1:1:10, StepRangeLen(1, 1, 10), LinRange(1, 10, 10), 0.1:0.2:0.9)
        s = OffsetArray(2:3, 2:3)
        @test axes(r[s]) == axes(s)
    end
end

@testset "logical indexing" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))

    @test A[A .> 2] == [3,4]
end

end