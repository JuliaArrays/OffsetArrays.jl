module OffsetArraysTest

using CatIndices: BidirectionalVector
using DelimitedFiles
using EllipsisNotation
using LinearAlgebra
using OffsetArraysCore
using OffsetArraysCore: no_offset_view, Origin, IdOffsetRange, IdentityUnitRange
using Test
using StaticArrays

include("customranges.jl")

no_offset_axes(x, d) = no_offset_view(axes(x, d))
no_offset_axes(x) = map(no_offset_view, axes(x))

function same_value(r1, r2)
    length(r1) == length(r2) || return false
    for (v1, v2) in zip(r1, r2)
        v1 == v2 || return false
    end
    return true
end

# Custom index types
struct ZeroBasedIndexing end
struct NewColon end
struct TupleOfRanges{N}
    x ::NTuple{N, UnitRange{Int}}
end

# used in testing the constructor
struct WeirdInteger{T} <: Integer
    x :: T
end
# assume that it doesn't behave as expected
Base.Int(a::WeirdInteger) = a

@testset "Constructors" begin
    @testset "Single-entry arrays in dims 0:5" begin
        for n = 0:5
            for z in (OffsetArray(ones(Int,ntuple(d->1,n)), ntuple(x->x-1,n)),
                    fill!(OffsetArray{Float64}(undef, ntuple(x->x:x, n)), 1),
                    fill!(OffsetArray{Float64}(undef, ntuple(x->x:x, n)...), 1),
                    fill!(OffsetArray{Float64,n}(undef, ntuple(x->x:x, n)), 1),
                    fill!(OffsetArray{Float64,n}(undef, ntuple(x->x:x, n)...), 1))
                @test length(LinearIndices(z)) == 1
                @test no_offset_axes(z) == ntuple(x->x:x, n)
                @test z[1] == 1
            end
        end
        a0 = reshape([3])
        a = OffsetArray(a0)
        @test axes(a) == ()
        @test ndims(a) == 0
        @test a[] == 3
        @test a === OffsetArray(a, ())
        @test_throws ArgumentError OffsetArray(a, 0)
        @test_throws ArgumentError OffsetArray(a0, 0)
    end

    @testset "OffsetVector" begin
        # initialization
        one_based_axes = Any[
            (Base.OneTo(4), ),
            (1:4, ),
            (big(1):big(4), ),
            (CartesianIndex(1):CartesianIndex(4), ),
            (IdentityUnitRange(1:4), ),
            (IdOffsetRange(1:4),),
            (IdOffsetRange(3:6, -2),)
        ]

        offset_axes = Any[
            (-1:2, ),
            (big(-1):big(2), ),
            (CartesianIndex(-1):CartesianIndex(2), ),
            (IdentityUnitRange(-1:2), ),
            (IdOffsetRange(-1:2),),
            (IdOffsetRange(3:6, -4),)
        ]

        offsets = size.(one_based_axes[1], 1)
        offsets_big = map(big, offsets)

        for inds in Any[offsets, offsets_big, one_based_axes...]
            # test indices API
            a = OffsetVector{Float64}(undef, inds)
            @test eltype(a) === Float64
            @test axes(a) === axes(OffsetVector{Float64}(undef, inds...)) === axes(OffsetArray{Float64, 1}(undef, inds)) === axes(OffsetArray{Float64}(undef, inds))
            @test axes(a) === (IdOffsetRange(Base.OneTo(4), 0), )
            @test a.offsets === (0, )
            @test axes(a.parent) == (Base.OneTo(4), )

            a = OffsetVector{Nothing}(nothing, inds)
            @test eltype(a) === Nothing
            @test axes(a) === axes(OffsetVector{Nothing}(nothing, inds...)) === axes(OffsetArray{Nothing, 1}(nothing, inds))
            @test axes(a) === (IdOffsetRange(Base.OneTo(4), 0), )

            a = OffsetVector{Missing}(missing, inds)
            @test eltype(a) === Missing
            @test axes(a) === axes(OffsetVector{Missing}(missing, inds...)) === axes(OffsetArray{Missing, 1}(missing, inds))
            @test axes(a) === (IdOffsetRange(Base.OneTo(4), 0), )
        end

        # nested OffsetVectors
        for inds in Any[offsets, offsets_big]
            a = OffsetVector{Float64}(undef, inds)
            b = OffsetVector(a, inds); b2 = OffsetVector(a, inds...);
            @test eltype(b) === eltype(b2) === Float64
            @test axes(b, 1) === axes(b2, 1) === IdOffsetRange(Base.OneTo(4), 4)
        end

        # offset indexing
        for inds in offset_axes
            # test offsets
            a = OffsetVector{Float64}(undef, inds)
            ax = (IdOffsetRange(Base.OneTo(4), -2), )
            @test a.offsets === (-2, )
            @test axes(a.parent) == (Base.OneTo(4), )
            @test axes(a) === ax
            a = OffsetVector{Nothing}(nothing, inds)
            @test axes(a) === ax
            a = OffsetVector{Missing}(missing, inds)
            @test axes(a) === ax

            for (T, t) in Any[(Nothing, nothing), (Missing, missing)]
                a = OffsetVector{Union{T, Vector{Int}}}(undef, inds)
                @test !isassigned(a, -1)
                @test eltype(a) === Union{T, Vector{Int}}
                @test axes(a) === ax

                a = OffsetVector{Union{T, Vector{Int}}}(t, inds)
                @test a[-1] === t
            end
        end
        @test_throws Union{ArgumentError, ErrorException} OffsetVector{Float64}(undef, -2) # only positive number works

        # convenient constructors
        a = rand(4)
        for inds in offset_axes
            oa1 = OffsetVector(a, inds...)
            oa2 = OffsetVector(a, inds)
            oa3 = OffsetArray(a, inds...)
            oa4 = OffsetArray(a, inds)
            @test oa1 === oa2 === oa3 === oa4
            @test axes(oa1) === (IdOffsetRange(Base.OneTo(4), -2), )
            @test parent(oa1) === a
            @test oa1.offsets === (-2, )
        end

        oa = OffsetArray(a, :)
        @test oa === OffsetArray(a, (:, )) === OffsetArray(a, axes(a)) === OffsetVector(a, :) === OffsetVector(a, axes(a))
        @test oa == a
        @test axes(oa) == axes(a)
        @test axes(oa) !== axes(a)

        # nested offset array
        a = rand(4)
        oa = OffsetArray(a, -1)
        for inds in Any[.-oa.offsets, one_based_axes...]
            ooa = OffsetArray(oa, inds)
            @test typeof(parent(ooa)) <: Vector
            @test ooa === OffsetArray(oa, inds...) === OffsetVector(oa, inds) === OffsetVector(oa, inds...)
            @test ooa == a
            @test axes(ooa) == axes(a)
            @test axes(ooa) !== axes(a)
        end

        # overflow bounds check
        v = rand(5)
        @test axes(OffsetVector(v, typemax(Int)-length(v))) == (IdOffsetRange(axes(v)[1], typemax(Int)-length(v)), )
        @test_throws OverflowError OffsetVector(v, typemax(Int)-length(v)+1)
        ao = OffsetArray(v, typemin(Int))
        ao2 = OffsetArray{Float64, 1, typeof(ao)}(ao, (-1, ))
        @test no_offset_axes(ao2, 1) == typemin(Int) .+ (0:length(v)-1)
        ao2 = OffsetArray(ao, (-1,))
        @test no_offset_axes(ao2, 1) == typemin(Int) .+ (0:length(v)-1)
        @test_throws OverflowError OffsetArray{Float64, 1, typeof(ao)}(ao, (-2, )) # inner Constructor
        @test_throws OverflowError OffsetArray(ao, (-2, )) # convenient constructor accumulate offsets
        @test_throws OverflowError OffsetVector(1:0, typemax(Int))
        @test_throws OverflowError OffsetVector(OffsetVector(1:0, 0), typemax(Int))
        @test_throws OverflowError OffsetArray(zeros(Int, typemax(Int):typemax(Int)), 2)
        @test_throws OverflowError OffsetArray(v, Origin(typemax(Int)))

        b = OffsetArray(OffsetArray(big(1):2, 1), typemax(Int)-1)
        @test no_offset_axes(b, 1) == big(typemax(Int)) .+ (1:2)

        @testset "OffsetRange" begin
            for r in Any[1:100, big(1):big(2)]
                a = OffsetVector(r, 4)
                @test first(r) in a
                @test !(last(r) + 1 in a)
            end

            @testset "BigInt axes" begin
                r = OffsetArray(1:big(2)^65, 4000)
                @test eltype(r) === BigInt
                @test no_offset_axes(r, 1) == (big(1):big(2)^65) .+ 4000
            end
        end

        # disallow OffsetVector(::Array{<:Any, N}, offsets) where N != 1
        @test_throws ArgumentError OffsetVector(zeros(2,2), (2, 2))
        @test_throws ArgumentError OffsetVector(zeros(2,2), 2, 2)
        @test_throws ArgumentError OffsetVector(zeros(2,2), (1:2, 1:2))
        @test_throws ArgumentError OffsetVector(zeros(2,2), 1:2, 1:2)
        @test_throws ArgumentError OffsetVector(zeros(), ())
        @test_throws ArgumentError OffsetVector(zeros())
        @test_throws ArgumentError OffsetVector(zeros(2,2), ())
        @test_throws ArgumentError OffsetVector(zeros(2,2))
        @test_throws ArgumentError OffsetVector(zeros(2,2), 2)
        @test_throws ArgumentError OffsetVector(zeros(2,2), (2,))
        @test_throws ArgumentError OffsetVector(zeros(2:3,2:3), 2, 3)
        @test_throws ArgumentError OffsetVector(zeros(2:3,2:3), (2, 4))
        @test_throws ArgumentError OffsetVector(zeros(2:3,2:3), ())
        @test_throws ArgumentError OffsetVector(zeros(2:3,2:3))

        # eltype of an OffsetArray should match that of the parent (issue #162)
        @test_throws TypeError OffsetVector{Float64,Vector{ComplexF64}}
        # ndim of an OffsetArray should match that of the parent
        @test_throws TypeError OffsetVector{Float64,Matrix{Float64}}
    end

    @testset "OffsetMatrix" begin
        # initialization

        one_based_axes = Any[
                (Base.OneTo(4), Base.OneTo(3)),
                (1:4, 1:3),
                (big(1):big(4), big(1):big(3)),
                (CartesianIndex(1, 1):CartesianIndex(4, 3), ),
                (CartesianIndex(1):CartesianIndex(4), CartesianIndex(1):CartesianIndex(3)),
                (CartesianIndex(1):CartesianIndex(4), 1:3),
                (IdentityUnitRange(1:4), IdentityUnitRange(1:3)),
                (IdOffsetRange(1:4), IdOffsetRange(1:3)),
                (IdOffsetRange(3:6, -2), IdOffsetRange(3:5, -2)),
                (IdOffsetRange(3:6, -2), IdentityUnitRange(1:3)),
                (IdOffsetRange(3:6, -2), 1:3),
        ]

        offset_axes = Any[
                (-1:2, 0:2),
                (big(-1):big(2), big(0):big(2)),
                (CartesianIndex(-1, 0):CartesianIndex(2, 2), ),
                (-1:2, CartesianIndex(0):CartesianIndex(2)),
                (CartesianIndex(-1):CartesianIndex(2), CartesianIndex(0):CartesianIndex(2)),
                (CartesianIndex(-1):CartesianIndex(2), 0:2),
                (IdentityUnitRange(-1:2), 0:2),
                (IdOffsetRange(-1:2), IdOffsetRange(0:2)),
                (IdOffsetRange(3:6, -4), IdOffsetRange(2:4, -2)),
                (IdOffsetRange(3:6, -4), IdentityUnitRange(0:2)),
                (IdOffsetRange(-1:2), 0:2),
        ]

        offsets = size.(one_based_axes[1], 1)
        offsets_big = map(big, offsets)

        for inds in Any[offsets, offsets_big, one_based_axes...]
            # test API
            a = OffsetMatrix{Float64}(undef, inds)
            ax = (IdOffsetRange(Base.OneTo(4), 0), IdOffsetRange(Base.OneTo(3), 0))
            @test eltype(a) === Float64
            @test axes(a) === axes(OffsetMatrix{Float64}(undef, inds...)) === axes(OffsetArray{Float64, 2}(undef, inds)) === axes(OffsetArray{Float64, 2}(undef, inds...)) === axes(OffsetArray{Float64}(undef, inds))
            @test axes(a) === ax
            @test a.offsets === (0, 0)
            @test axes(a.parent) == (Base.OneTo(4), Base.OneTo(3))

            a = OffsetMatrix{Nothing}(nothing, inds)
            @test eltype(a) === Nothing
            @test axes(a) === axes(OffsetMatrix{Nothing}(nothing, inds...)) === axes(OffsetArray{Nothing, 2}(nothing, inds)) === axes(OffsetArray{Nothing, 2}(nothing, inds...))
            @test axes(a) === ax

            a = OffsetMatrix{Missing}(missing, inds)
            @test eltype(a) === Missing
            @test axes(a) === axes(OffsetMatrix{Missing}(missing, inds...)) === axes(OffsetArray{Missing, 2}(missing, inds)) === axes(OffsetArray{Missing, 2}(missing, inds...))
            @test axes(a) === ax
        end
        @test_throws Union{ArgumentError, ErrorException} OffsetMatrix{Float64}(undef, 2, -2) # only positive numbers works

        # nested OffsetMatrices
        for inds in Any[offsets, offsets_big]
            a = OffsetMatrix{Float64}(undef, inds)
            b = OffsetMatrix(a, inds); b2 = OffsetMatrix(a, inds...);
            @test eltype(b) === eltype(b2) === Float64
            @test axes(b, 1) === axes(b2, 1) === IdOffsetRange(Base.OneTo(4), 4)
            @test axes(b, 2) === axes(b2, 2) === IdOffsetRange(Base.OneTo(3), 3)
        end

        for inds in offset_axes
            # test offsets
            a = OffsetMatrix{Float64}(undef, inds)
            ax = (IdOffsetRange(Base.OneTo(4), -2), IdOffsetRange(Base.OneTo(3), -1))
            @test a.offsets === (-2, -1)
            @test axes(a.parent) == (Base.OneTo(4), Base.OneTo(3))
            @test axes(a) === ax
            a = OffsetMatrix{Nothing}(nothing, inds)
            @test axes(a) === ax
            a = OffsetMatrix{Missing}(missing, inds)
            @test axes(a) === ax

            for (T, t) in Any[(Nothing, nothing), (Missing, missing)]
                a = OffsetMatrix{Union{T, Vector{Int}}}(undef, inds)
                @test !isassigned(a, -1, 0)
                @test eltype(a) === Union{T, Vector{Int}}
                @test axes(a) === ax

                a = OffsetMatrix{Union{T, Vector{Int}}}(t, inds)
                @test a[-1, 0] === t
            end
        end

        # convenient constructors
        a = rand(4, 3)
        for inds in offset_axes
            ax = (IdOffsetRange(Base.OneTo(4), -2), IdOffsetRange(Base.OneTo(3), -1))
            oa1 = OffsetMatrix(a, inds...)
            oa2 = OffsetMatrix(a, inds)
            oa3 = OffsetArray(a, inds...)
            oa4 = OffsetArray(a, inds)
            @test oa1 === oa2 === oa3 === oa4
            @test axes(oa1) === ax
            @test parent(oa1) === a
            @test oa1.offsets === (-2, -1)
        end
        oa = OffsetArray(a, :, axes(a, 2))
        @test oa === OffsetArray(a, (axes(oa, 1), :)) === OffsetArray(a, axes(a)) === OffsetMatrix(a, (axes(oa, 1), :)) === OffsetMatrix(a, axes(a))
        @test oa == a
        @test axes(oa) == axes(a)
        @test axes(oa) !== axes(a)

        oa = OffsetMatrix(a, :, 2:4)
        @test oa === OffsetMatrix(a, axes(a, 1), 2:4) === OffsetMatrix(a, (axes(oa, 1), 2:4))

        # nested offset array
        a = rand(4, 3)
        oa = OffsetArray(a, -1, -2)
        for inds in Any[.-oa.offsets, one_based_axes...]
            ooa = OffsetArray(oa, inds)
            @test ooa === OffsetArray(oa, inds...) === OffsetMatrix(oa, inds) === OffsetMatrix(oa, inds...)
            @test typeof(parent(ooa)) <: Matrix
            @test ooa == a
            @test axes(ooa) == axes(a)
            @test axes(ooa) !== axes(a)
        end

        # overflow bounds check
        a = rand(4, 3)
        @test axes(OffsetMatrix(a, typemax(Int)-size(a, 1), 0)) == (IdOffsetRange(axes(a)[1], typemax(Int)-size(a, 1)), axes(a, 2))
        @test_throws OverflowError OffsetMatrix(a, typemax(Int)-size(a,1)+1, 0)
        @test_throws OverflowError OffsetMatrix(a, 0, typemax(Int)-size(a, 2)+1)

        # disallow OffsetMatrix(::Array{<:Any, N}, offsets) where N != 2
        @test_throws ArgumentError OffsetMatrix(zeros(2), (2,))
        @test_throws ArgumentError OffsetMatrix(zeros(2), 2)
        @test_throws ArgumentError OffsetMatrix(zeros(2), (1:2,))
        @test_throws ArgumentError OffsetMatrix(zeros(2), 1:2)
        @test_throws ArgumentError OffsetMatrix(zeros(), ())
        @test_throws ArgumentError OffsetMatrix(zeros())
        @test_throws ArgumentError OffsetMatrix(zeros(2), ())
        @test_throws ArgumentError OffsetMatrix(zeros(2))
        @test_throws ArgumentError OffsetMatrix(zeros(2), (1, 2))
        @test_throws ArgumentError OffsetMatrix(zeros(2), 1, 2)
        @test_throws ArgumentError OffsetMatrix(zeros(2:3), (2,))
        @test_throws ArgumentError OffsetMatrix(zeros(2:3), 2)
        @test_throws ArgumentError OffsetMatrix(zeros(2:3, 1:2, 1:2), (2,0,0))
        @test_throws ArgumentError OffsetMatrix(zeros(2:3, 1:2, 1:2), 2,0,0)
        @test_throws ArgumentError OffsetMatrix(zeros(2:3, 1:2, 1:2), ())
        @test_throws ArgumentError OffsetMatrix(zeros(2:3, 1:2, 1:2))

        # eltype of an OffsetArray should match that of the parent (issue #162)
        @test_throws TypeError OffsetMatrix{Float64,Matrix{ComplexF64}}
        # ndim of an OffsetArray should match that of the parent
        @test_throws TypeError OffsetMatrix{Float64,Vector{Float64}}
    end

    # no need to duplicate the 2D case here,
    # only add some special test cases
    @testset "OffsetArray" begin
        a = rand(2, 2, 2)
        oa = OffsetArray(a, 0:1, 3:4, 2:3)
        @test OffsetArray(a, CartesianIndices(axes(oa))) == oa
        @test no_offset_axes(OffsetArray(a, :, CartesianIndices((3:4, 2:3)))) == (1:2, 3:4, 2:3)
        @test no_offset_axes(OffsetArray(a, 10:11, CartesianIndices((3:4, 2:3)) )) == (10:11, 3:4, 2:3)
        @test no_offset_axes(OffsetArray(a, CartesianIndices((3:4, 2:3)), :)) == (3:4, 2:3, 1:2)
        @test no_offset_axes(OffsetArray(a, CartesianIndices((3:4, 2:3)), 10:11)) == (3:4, 2:3, 10:11)
        @test no_offset_axes(OffsetArray(a, :, :, CartesianIndices((3:4,)) )) == (1:2, 1:2, 3:4)
        @test no_offset_axes(OffsetArray(a, 10:11, :, CartesianIndices((3:4,)) )) == (10:11, 1:2, 3:4)
        @test no_offset_axes(OffsetArray(a, 10:11, 2:3, CartesianIndices((3:4,)) )) == (10:11, 2:3, 3:4)

        # ignore empty CartesianIndices
        @test OffsetArray(a, CartesianIndices(()), 0:1, :, 2:3) == OffsetArray(a, 0:1, :, 2:3)
        @test OffsetArray(a, 0:1, CartesianIndices(()), :, 2:3) == OffsetArray(a, 0:1, :, 2:3)
        @test OffsetArray(a, 0:1, :,  CartesianIndices(()), 2:3) == OffsetArray(a, 0:1, :, 2:3)
        @test OffsetArray(a, 0:1, :, 2:3, CartesianIndices(())) == OffsetArray(a, 0:1, :, 2:3)

        # nested OffsetArrays
        for offsets in [(1,1,1), big.((1,1,1))]
            ob = OffsetArray(oa, offsets); ob2 = OffsetArray(oa, offsets...);
            @test eltype(ob) === eltype(ob2) === Float64
            @test axes(ob, 1) === axes(ob2, 1) === IdOffsetRange(Base.OneTo(2), 0)
            @test axes(ob, 2) === axes(ob2, 2) === IdOffsetRange(Base.OneTo(2), 3)
            @test axes(ob, 3) === axes(ob2, 3) === IdOffsetRange(Base.OneTo(2), 2)
        end

        indices = (-1:1, -7:7, -1:2, -5:5, -1:1, -3:3, -2:2, -1:1)
        y = OffsetArray{Float64}(undef, indices...);
        @test axes(y) === axes(OffsetArray{Float64}(undef, indices))
        @test axes(y) === axes(OffsetArray{Float64, length(indices)}(undef, indices...))
        @test no_offset_axes(y) == (-1:1, -7:7, -1:2, -5:5, -1:1, -3:3, -2:2, -1:1)
        @test eltype(y) === Float64

        @test_throws ArgumentError OffsetArray{Float64, 2}(undef, indices)
        @test_throws ArgumentError OffsetArray(y, indices[1:2])

        @test ndims(OffsetArray(zeros(), ())) == 0
        @test Base.axes1(OffsetArray(zeros(), ())) === IdOffsetRange(Base.OneTo(1))

        # eltype of an OffsetArray should match that of the parent (issue #162)
        @test_throws TypeError OffsetArray{Float64,2,Matrix{ComplexF64}}
        # ndim of an OffsetArray should match that of the parent
        @test_throws TypeError OffsetArray{Float64,3,Matrix{Float64}}

        # should throw a TypeError if the offsets can not be converted to Ints
        @test_throws TypeError OffsetVector{Int,Vector{Int}}(zeros(Int,2), (WeirdInteger(1),))
    end

    @testset "custom range types" begin
        @testset "EllipsisNotation" begin
            @testset "Vector" begin
                v = rand(5)
                @test axes(OffsetArray(v, ..)) == axes(v)
                @test OffsetArray(v, ..) == OffsetArray(v, :)
                @test axes(OffsetVector(v, ..)) == axes(v)
                @test OffsetVector(v, ..) == OffsetVector(v, :)

                @test no_offset_axes(OffsetArray(v, .., 2:6)) == (2:6, )
                @test OffsetArray(v, .., 2:6) == OffsetArray(v, 2:6)
                @test no_offset_axes(OffsetVector(v, .., 2:6)) == (2:6, )
                @test OffsetVector(v, .., 2:6) == OffsetVector(v, 2:6)
            end
            @testset "Matrix" begin
                m = rand(2, 2)
                @test axes(OffsetArray(m, ..)) == axes(m)
                @test OffsetArray(m, ..) == OffsetArray(m, :, :)
                @test axes(OffsetMatrix(m, ..)) == axes(m)
                @test OffsetMatrix(m, ..) == OffsetMatrix(m, :, :)

                @test no_offset_axes(OffsetArray(m, .., 2:3)) == (axes(m, 1), 2:3)
                @test OffsetArray(m, .., 2:3) == OffsetArray(m, :, 2:3)
                @test no_offset_axes(OffsetMatrix(m, .., 2:3)) == (axes(m, 1), 2:3)
                @test OffsetMatrix(m, .., 2:3) == OffsetMatrix(m, :, 2:3)

                @test no_offset_axes(OffsetArray(m, .., 2:3, 3:4)) == (2:3, 3:4)
                @test OffsetArray(m, .., 2:3, 3:4) == OffsetArray(m, 2:3, 3:4)
                @test no_offset_axes(OffsetMatrix(m, .., 2:3, 3:4)) == (2:3, 3:4)
                @test OffsetMatrix(m, .., 2:3, 3:4) == OffsetMatrix(m, 2:3, 3:4)
            end
            @testset "3D Array" begin
                a = rand(2, 2, 2)
                @test axes(OffsetArray(a, ..)) == axes(a)
                @test OffsetArray(a, ..) == OffsetArray(a, :, :, :)

                @test no_offset_axes(OffsetArray(a, .., 2:3)) == (axes(a)[1:2]..., 2:3)
                @test OffsetArray(a, .., 2:3) == OffsetArray(a, :, :, 2:3)

                @test no_offset_axes(OffsetArray(a, .., 2:3, 3:4)) == (axes(a, 1), 2:3, 3:4)
                @test OffsetArray(a, .., 2:3, 3:4) == OffsetArray(a, :, 2:3, 3:4)

                @test no_offset_axes(OffsetArray(a, 2:3, .., 3:4)) == (2:3, axes(a, 2), 3:4)
                @test OffsetArray(a, 2:3, .., 3:4) == OffsetArray(a, 2:3, :, 3:4)

                @test no_offset_axes(OffsetArray(a, .., 4:5, 2:3, 3:4)) == (4:5, 2:3, 3:4)
                @test OffsetArray(a, .., 4:5, 2:3, 3:4) == OffsetArray(a, 4:5, 2:3, 3:4)
            end
        end
        @testset "ZeroBasedIndexing" begin
            Base.to_indices(A, inds, ::Tuple{ZeroBasedIndexing}) = map(x -> 0:length(x) - 1, inds)

            a = zeros(3,3)
            oa = OffsetArray(a, ZeroBasedIndexing())
            @test no_offset_axes(oa) == (0:2, 0:2)
        end
        @testset "TupleOfRanges" begin
            Base.to_indices(A, inds, t::Tuple{TupleOfRanges{N}}) where {N} = t
            OffsetArraysCore.AxisConversionStyle(::Type{TupleOfRanges{N}}) where {N} =
                OffsetArraysCore.TupleOfRanges()

            Base.convert(::Type{Tuple{Vararg{AbstractUnitRange{Int}}}}, t::TupleOfRanges) = t.x

            a = zeros(3,3)
            inds = TupleOfRanges((3:5, 2:4))
            oa = OffsetArray(a, inds)
            @test no_offset_axes(oa) == inds.x
        end
        @testset "NewColon" begin
            Base.to_indices(A, inds, t::Tuple{NewColon,Vararg{Any}}) =
                (_uncolon(inds, t), to_indices(A, Base.tail(inds), Base.tail(t))...)

            _uncolon(inds::Tuple{}, I::Tuple{NewColon, Vararg{Any}}) = OneTo(1)
            _uncolon(inds::Tuple, I::Tuple{NewColon, Vararg{Any}}) = inds[1]

            a = zeros(3, 3)
            oa = OffsetArray(a, (NewColon(), 2:4))
            @test no_offset_axes(oa) == (axes(a,1), 2:4)
        end
    end

    @testset "Offset range construction" begin
        r = -2:5
        for AT in Any[OffsetArray, OffsetVector]
            y = AT(r, r)
            @test no_offset_axes(y) == (r,)
            @test step(y) == step(r)
            y = AT(r, (r,))
            @test no_offset_axes(y) == (r,)
            y = AT(r, CartesianIndices((r, )))
            @test no_offset_axes(y) == (r, )
        end
    end

    @testset "size/length" begin
        for p in Any[SA[1,2,3,4], 1:4, [1:4;]]
            for A in Any[OffsetArray(p, 4),
                    OffsetArray(reshape(p, 2, 2), 3, 4),
                    OffsetArray(reshape(p, 2, 1, 2), 3, 0, 4),
                    OffsetArray(reshape(p, Val(1)), 2)]
                @test size(A) == size(parent(A))
                @test length(A) == length(parent(A))
            end
        end
    end
end

@testset "Axes supplied to constructor correspond to final result" begin
    # Ref https://github.com/JuliaArrays/OffsetArrays.jl/pull/65#issuecomment-457181268
    B = BidirectionalVector([1, 2, 3], -2)
    A = OffsetArray(B, -1:1)
    @test no_offset_axes(A) == (-1:1,)
end

@testset "view" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))

    S = view(A, :, 3)
    @test S == OffsetArray([1,2], (A.offsets[1],))
    @test S[0] == 1
    @test S[1] == 2
    @test_throws BoundsError S[2]
    @test axes(S) == (IdentityUnitRange(0:1),)
    S = view(A, 0, :)
    @test S == OffsetArray([1,3], (A.offsets[2],))
    @test S[3] == 1
    @test S[4] == 3
    @test_throws BoundsError S[1]
    @test axes(S) == (IdentityUnitRange(3:4),)
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
    @test axes(S) == IdentityUnitRange.((0:1, 3:4))
    S = view(A, axes(A)...)
    @test S == A
    @test S[0,3] == S[1] == 1
    @test S[1,3] == S[2] == 2
    @test S[0,4] == S[3] == 3
    @test S[1,4] == S[4] == 4
    @test_throws BoundsError S[1,1]
    @test no_offset_axes(S) == (0:1, 3:4)
    # issue 100
    S = view(A, axes(A, 1), 3)
    @test S == A[:, 3]
    @test S[0] == 1
    @test S[1] == 2
    @test_throws BoundsError S[length(S)]
    @test no_offset_axes(S) == (0:1, )
    # issue 100
    S = view(A, 1, axes(A, 2))
    @test S == A[1, :]
    @test S[3] == 2
    @test S[4] == 4
    @test_throws BoundsError S[1]
    @test no_offset_axes(S) == (3:4, )

    # issue 133
    r = IdOffsetRange(1:2, -1)
    v1 = view(A, r, 3)
    @test v1[0] == 1
    @test v1[1] == 2
    @test axes(v1, 1) == axes(r, 1)
    v2 = view(A, UnitRange(r), 3)
    for (indflat, indoffset) in enumerate(r)
        @test v1[indoffset] == v2[indflat]
    end

    # issue 133
    r = IdOffsetRange(1:2, 2)
    v1 = view(A, 1, r)
    @test v1[3] == 2
    @test v1[4] == 4
    @test axes(v1, 1) == axes(r, 1)
    v2 = view(A, 1, UnitRange(r))
    for (indflat, indoffset) in enumerate(r)
        @test v1[indoffset] == v2[indflat]
    end

    # issue 133
    a12 = zeros(3:8, 3:4)
    r = IdOffsetRange(Base.OneTo(3), 5)
    a12[r, 4] .= 3
    @test all(a12[r, 4] .== 3)
    @test all(a12[UnitRange(r), 4] .== 3)

    A0 = collect(reshape(1:24, 2, 3, 4))
    A = OffsetArray(A0, (-1,2,1))
    S = view(A, axes(A, 1), 3:4, axes(A, 3))
    @test S == A[:, 3:4, :]
    @test S[0, 1, 2] == A[0, 3, 2]
    @test S[0, 2, 2] == A[0, 4, 2]
    @test S[1, 1, 2] == A[1, 3, 2]
    @test no_offset_axes(S) == (0:1, Base.OneTo(2), 2:5)

    # issue #186
    a = reshape(1:12, 3, 4)
    r = IdOffsetRange(3:4)
    av = view(a, :, r)
    @test av == a[:, 3:4]
    @test axes(av) == (axes(a,1), axes(r,1))
    r = IdOffsetRange(1:2,2)
    av = view(a, :, r)
    @test no_offset_view(av) == a[:, 3:4]
    @test axes(av) == (axes(a,1), axes(r,1))
    r = IdOffsetRange(2:3)
    av1d = view(a, r, 3)
    @test av1d == a[2:3, 3]
    @test axes(av1d) == (axes(r,1),)
    r = IdOffsetRange(Base.OneTo(2), 1)
    av1d = view(a, r, 3)
    @test no_offset_view(av1d) == a[2:3, 3]
    @test axes(av1d) == (axes(r,1),)

    # fix IdOffsetRange(::IdOffsetRange, offset) nesting from #178
    b = 1:20
    bov = OffsetArray(view(b, 3:4), 3:4)
    c = @view b[bov]
    @test same_value(c, 3:4)
    @test no_offset_axes(c,1) == 3:4
    d = OffsetArray(c, 1:2)
    @test same_value(d, c)
    @test axes(d,1) == 1:2

    # Issue 128
    a = OffsetArray(1:3, 0:2);
    b = @view a[0]
    @test b[] == b[1] == 1

    a = reshape(1:16, 4, 4);
    for ax1 in (:, axes(a,1), UnitRange(axes(a,1))),
        ax2 in (:, axes(a,2), UnitRange(axes(a,2)))
            av = @view a[ax1, ax2]
            av_nooffset = no_offset_view(av)
            @test axes(av_nooffset) === axes(av)
    end
end


@testset "iteration" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))

    let a
        for (a,d) in zip(A, A0)
            @test a == d
        end
    end

    v = ones(10)
    for r in Any[1:1:10, 1:10], s in Any[r, collect(r)]
        so = OffsetArray(s)
        @test Float64[v[i] for i in s] == Float64[v[i] for i in so]
    end
end

@testset "show/summary" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))
    S = OffsetArray(view(A0, 1:2, 1:2), (-1,2))

    @test sprint(show, A) == "[1 3; 2 4]"
    @test sprint(show, S) == "[1 3; 2 4]"
    strs = split(strip(sprint(show, MIME("text/plain"), A)), '\n')
    @test strs[2] == " 1  3"
    @test strs[3] == " 2  4"
    v = OffsetArray(rand(3), (-2,))
    @test sprint(show, v) == sprint(show, parent(v))
    function cmp_showf(printfunc, A)
        str1 = sprint(printfunc, A, context=(:limit=>true, :compact=>true))
        str2 = sprint(printfunc, parent(A), context=(:limit=>true, :compact=>true))
        @test str1 == str2
    end
    cmp_showf(Base.print_matrix, OffsetArray(rand(5,5), (10,-9)))       # rows&cols fit
    cmp_showf(Base.print_matrix, OffsetArray(rand(10^3,5), (10,-9)))    # columns fit
    cmp_showf(Base.print_matrix, OffsetArray(rand(5,10^3), (10,-9)))    # rows fit
    cmp_showf(Base.print_matrix, OffsetArray(rand(10^3,10^3), (10,-9))) # neither fits

    a = OffsetArray([1 2; 3 4], -1:0, 5:6)
    shownsz = VERSION >= v"1.2.0-DEV.229" ? Base.dims2string(size(a))*' ' : ""
    @test summary(a) == "$(shownsz)OffsetArray(::$(typeof(parent(a))), -1:0, 5:6) with eltype $(Int) with indices -1:0×5:6"
    shownsz = VERSION >= v"1.2.0-DEV.229" ? Base.dims2string(size(view(a, :, 5)))*' ' : ""
    @test summary(view(a, :, 5)) == "$(shownsz)view(OffsetArray(::$(typeof(parent(a))), -1:0, 5:6), :, 5) with eltype $(Int) with indices -1:0"
    a = OffsetArray(reshape([1]))
    @test summary(a) == "0-dimensional OffsetArray(::$(typeof(parent(a)))) with eltype $(Int)"

    a = OffsetArray([1 2; 3 4], -1:0, 5:6)
    @test occursin("IdOffsetRange(values=-1:0, indices=-1:0)", sprint(show, axes(a, 1)))
    @test occursin("IdOffsetRange(values=5:6, indices=5:6)", sprint(show, axes(a, 2)))
    rrtable = IdOffsetRange(values=7:9, indices=-1:1)
    rrted = eval(Meta.parse(string(rrtable)))
    @test pairs(rrtable) == pairs(rrted)

    @test Base.inds2string(axes(a)) == Base.inds2string(map(UnitRange, axes(a)))
    @test Base.inds2string((IdOffsetRange(3:4),)) == "3:4"
    @test Base.inds2string((IdentityUnitRange(IdOffsetRange(3:4)),)) == "3:4"
    # check that the following doesn't throw
    @test Base.inds2string(()) isa Any

    @test sprint(show, OffsetArray(3:5, 0:2)) == "3:5 with indices 0:2"

    @test sprint(show, MIME"text/plain"(), OffsetArray(3:5, 0:2)) == "3:5 with indices 0:2"

    # issue #198
    for r in Any[axes(OffsetVector(1:10, -5), 1), 1:1:2, 1.0:1.0:2.0, 1:-1:-5]
        a = OffsetVector(r, 5)
        @test sprint(show, a) == "$r with indices $(UnitRange(axes(a,1)))"
    end

    d = Diagonal([1,2,3])
    s1 = sprint(Base.print_array, d)
    od = OffsetArray(d, -1:1, 3:5)
    s2 = sprint(Base.print_array, od)
    @test s1 == s2

    @test Base.replace_in_print_matrix(od, -1, 3, " ") == Base.replace_in_print_matrix(d, 1, 1, " ")
    @test Base.replace_in_print_matrix(od, -1, 4, " ") == Base.replace_in_print_matrix(d, 1, 2, " ")

    v = rand(3)
    ov = OffsetArray(v, (-2,))
    @test Base.replace_in_print_matrix(ov, -1, 1, " ") == Base.replace_in_print_matrix(v, 1, 1, " ")

    # Avoid returning the value of toplevel if it is false
    # showarg should only print values, it shouldn't return anything
    io = IOBuffer()
    @test Base.showarg(io, a, false) === nothing
    # check the other case too for good measure
    @test Base.showarg(io, a, true) === nothing
end

@testset "readdlm/writedlm" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))

    io = IOBuffer()
    writedlm(io, A)
    seek(io, 0)
    @test readdlm(io, eltype(A)) == parent(A)
end

@testset "similar" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))

    B = similar(A, Float32)
    @test isa(B, OffsetArray{Float32,2})
    @test axes(B) === axes(A)
    B = similar(A, (3,4))
    @test isa(B, Array{Int,2})
    @test size(B) == (3,4)
    @test axes(B) === (Base.OneTo(3), Base.OneTo(4))
    B = similar(A, (-3:3,1:4))
    @test isa(B, OffsetArray{Int,2})
    @test axes(B) == IdentityUnitRange.((-3:3, 1:4))

    idunitrange(x::AbstractUnitRange) = IdOffsetRange(values=x, indices=x)
    idunitrange(x::Integer) = x

    B = similar(parent(A), idunitrange.((-3:3,1:4)))
    @test isa(B, OffsetArray{Int,2})
    @test axes(B) == IdentityUnitRange.((-3:3, 1:4))
    @test isa([x for x in [1,2,3]], Vector{Int})
    @test similar(Array{Int}, idunitrange.((0:0, 0:0))) isa OffsetArray{Int, 2}
    @test similar(Array{Int}, (1, 1)) isa Matrix{Int}
    @test similar(Array{Int}, (Base.OneTo(1), Base.OneTo(1))) isa Matrix{Int}

    # check with an unseen axis type
    A = similar(ones(1), Int, ZeroBasedUnitRange(3:4), 4:5)
    @test eltype(A) === Int
    @test no_offset_axes(A) == (3:4, 4:5)
    A = similar(ones(1), Int, IdOffsetRange(ZeroBasedUnitRange(3:4), 2), 4:5)
    @test eltype(A) === Int
    @test no_offset_axes(A) == (5:6, 4:5)
    A = similar(ones(1), Int, IdOffsetRange(ZeroBasedUnitRange(3:4), 2), 4)
    @test eltype(A) === Int
    @test no_offset_axes(A) == (5:6, 1:4)

    @testset "similar with OffsetArray type (issue #263)" begin
        for i in Any[[1,2,3], 1:3, reshape(1:4, 2, 2)]
            k = OffsetArray(i, map(x -> -2, size(i)))
            j = similar(typeof(k), axes(k))
            @test axes(j) == axes(k)
            @test eltype(j) == eltype(k)
            j = similar(typeof(k), size(k))
            @test eltype(j) == eltype(k)
            @test size(j) == size(k)
            @test all(==(1), first.(axes(j)))
        end
    end
end

@testset "copyto!" begin
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
    if VERSION < v"1.5-"
        @test_throws BoundsError copyto!(a, b)
        fill!(a, -1)
        copyto!(a, bo)
        @test a[-3] == -1
        @test a[-2] == 1
        @test a[-1] == 2
    else
        # the behavior of copyto! is corrected as the documentation says "first n element"
        # https://github.com/JuliaLang/julia/pull/34049
        fill!(a, -1)
        copyto!(a, bo)
        @test a[-3] == 1
        @test a[-2] == 2
        @test a[-1] == -1
    end
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
end

@testset "map" begin
    am = OffsetArray{Int}(undef, (1:1, 7:9))  # for testing linear indexing
    fill!(am, -1)
    copyto!(am, 1:2)

    dest = similar(am)
    map!(+, dest, am, am)
    @test dest[1,7] == 2
    @test dest[1,8] == 4
    @test dest[1,9] == -2

    @testset "eltype conversion" begin
        a = OffsetArray(1:2, 1)
        b = map(BigInt, a)
        @test eltype(b) == BigInt
        @test b == a
        @test parent(b) isa AbstractRange

        for ri in Any[2:3, Base.OneTo(2)]
            for r in Any[IdOffsetRange(ri), IdOffsetRange(ri, 1), OffsetArray(ri), OffsetArray(ri, 2)]
                for T in [Int8, Int16, Int32, Int64, Int128, BigInt, Float32, Float64, BigFloat]
                    r2 = map(T, r)
                    @test eltype(r2) == T
                    @test axes(r2) == axes(r)
                    @test all(((x,y),) -> isequal(x,y), zip(r, r2))
                end
            end
        end

        @testset "Bool" begin
            for ri in Any[0:0, 0:1, 1:0, 1:1, Base.OneTo(0), Base.OneTo(1)]
                for r = Any[IdOffsetRange(ri), IdOffsetRange(ri .- 1, 1), OffsetVector(ri)]
                    r2 = map(Bool, r)
                    @test eltype(r2) == Bool
                    @test axes(r2) == axes(r)
                    @test all(((x,y),) -> isequal(x,y), zip(r, r2))
                end
            end
        end
    end
end

@testset "reductions" begin
    A = OffsetArray(rand(Int,4,4), (-3,5))
    @test maximum(A) == maximum(parent(A))
    @test minimum(A) == minimum(parent(A))
    @test extrema(A) == extrema(parent(A))
    @test sum(A) == sum(parent(A))
    @test sum(A, dims=1) == OffsetArray(sum(parent(A), dims=1), A.offsets)
    @test sum(A, dims=2) == OffsetArray(sum(parent(A), dims=2), A.offsets)
    @test sum(A, dims=(1,2)) == OffsetArray(sum(parent(A), dims=(1,2)), A.offsets)
    @test sum(view(OffsetArray(reshape(1:27, 3, 3, 3), 0, 0, 0), :, :, 1:2), dims=(2,3)) == reshape([51,57,63], 3, 1, 1)
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

    amin, amax = extrema(parent(A))
    @test clamp.(A, (amax+amin)/2, amax) == OffsetArray(clamp.(parent(A), (amax+amin)/2, amax), axes(A))

    @testset "mapreduce for OffsetRange" begin
        rangelist = Any[
            # AbstractUnitRanges
            5:100, UnitRange(5.0, 20.0), false:true,
            IdOffsetRange(4:5),
            IdOffsetRange(1:10, 4),
            # AbstractRanges
            2:4:14, 1.5:1.0:10.5,
            ]

        for r in rangelist

            a = OffsetVector(r, 2);
            @test mapreduce(identity, +, a) == mapreduce(identity, +, r)
            @test mapreduce(x -> x^2, (x,y) -> x, a) == mapreduce(x -> x^2, (x,y) -> x, r)

            b = mapreduce(identity, +, a, dims = 1)
            br = mapreduce(identity, +, r, dims = 1)
            @test no_offset_view(b) == no_offset_view(br)
            @test no_offset_axes(b, 1) == first(axes(a,1)):first(axes(a,1))

            @test mapreduce(identity, +, a, init = 3) == mapreduce(identity, +, r, init = 3)
            if VERSION >= v"1.2"
                @test mapreduce((x,y) -> x*y, +, a, a) == mapreduce((x,y) -> x*y, +, r, r)
                @test mapreduce((x,y) -> x*y, +, a, a, init = 10) == mapreduce((x,y) -> x*y, +, r, r, init = 10)
            end

            for f in [sum, minimum, maximum]
                @test f(a) == f(r)

                b = f(a, dims = 1);
                br = f(r, dims = 1)
                @test no_offset_view(b) == no_offset_view(br)
                @test no_offset_axes(b, 1) == first(axes(a,1)):first(axes(a,1))

                b = f(a, dims = 2);
                br = f(r, dims = 2)
                @test no_offset_view(b) == no_offset_view(br)
                @test axes(b, 1) == axes(a,1)
            end

            @test extrema(a) == extrema(r)
        end
    end
end

@testset "Collections" begin
    A = OffsetArray(rand(4,4), (-3,5))

    @test unique(A, dims=1) == OffsetArray(parent(A), 0, first(axes(A, 2)) - 1)
    @test unique(A, dims=2) == OffsetArray(parent(A), first(axes(A, 1)) - 1, 0)
    v = OffsetArray(rand(8), (-2,))
    @test sort(v) == OffsetArray(sort(parent(v)), v.offsets)
    @test sortslices(A; dims=1) == OffsetArray(sortslices(parent(A); dims=1), A.offsets)
    @test sortslices(A; dims=2) == OffsetArray(sortslices(parent(A); dims=2), A.offsets)
    @test sort(A, dims = 1) == OffsetArray(sort(parent(A), dims = 1), A.offsets)
    @test sort(A, dims = 2) == OffsetArray(sort(parent(A), dims = 2), A.offsets)

    @test mapslices(sort, A, dims = 1) == OffsetArray(mapslices(sort, parent(A), dims = 1), A.offsets)
    @test mapslices(sort, A, dims = 2) == OffsetArray(mapslices(sort, parent(A), dims = 2), A.offsets)
end

@testset "rot/reverse" begin
    A = OffsetArray(rand(4,4), (-3,5))

    @test rotl90(A) == OffsetArray(rotl90(parent(A)), A.offsets[[2,1]])
    @test rotr90(A) == OffsetArray(rotr90(parent(A)), A.offsets[[2,1]])
    @test reverse(A, dims = 1) == OffsetArray(reverse(parent(A), dims = 1), A.offsets)
    @test reverse(A, dims = 2) == OffsetArray(reverse(parent(A), dims = 2), A.offsets)
end

@testset "broadcasting" begin
    A = OffsetArray(rand(4,4), (-3,5))

    @test A.+1 == OffsetArray(parent(A).+1, A.offsets)
    @test 2*A == OffsetArray(2*parent(A), A.offsets)
    @test A+A == OffsetArray(parent(A)+parent(A), A.offsets)
    @test A.*A == OffsetArray(parent(A).*parent(A), A.offsets)
end

@testset "Resizing OffsetVectors" begin
    local a = OffsetVector(rand(5),-3)
    axes(a,1) == -2:2
    length(a) == 5
    resize!(a,3)
    length(a) == 3
    axes(a,1) == -2:0
    @test_throws ArgumentError resize!(a,-3)
end

####
#### type defined for testing no_offset_view
####

struct NegativeArray{T,N,S <: AbstractArray{T,N}} <: AbstractArray{T,N}
    parent::S
end

# Note: this defines the axes-of-the-axes to be OneTo.
# In general this isn't recommended, because
#    positionof(A, i, j, ...) == map(getindex, axes(A), (i, j, ...))
# is quite desirable, and this requires that the axes be "identity" ranges, i.e.,
# `r[i] == i`.
# Nevertheless it's useful to test this on a "broken" implementation
# to make sure we still get the right answer.
Base.axes(A::NegativeArray) = map(n -> (-n):(-1), size(A.parent))

Base.size(A::NegativeArray) = size(A.parent)

function Base.getindex(A::NegativeArray{T,N}, I::Vararg{Int,N}) where {T,N}
    getindex(A.parent, (I .+ size(A.parent) .+ 1)...)
end

struct PointlessWrapper{T,N, A <: AbstractArray{T,N}} <: AbstractArray{T,N}
    parent :: A
end
Base.parent(x::PointlessWrapper) = x.parent
Base.size(x::PointlessWrapper) = size(parent(x))
Base.axes(x::PointlessWrapper) = axes(parent(x))
Base.getindex(x::PointlessWrapper, i...) = x.parent[i...]

@testset "no offset view" begin
    # OffsetArray fallback
    A = randn(3, 3)
    @inferred no_offset_view(A)
    O1 = OffsetArray(A, -1:1, 0:2)
    O2 = OffsetArray(O1, -2:0, -3:(-1))
    @test no_offset_view(O2) ≡ A
    @inferred no_offset_view(O1)
    @inferred no_offset_view(O2)

    P = PointlessWrapper(A)
    @test @inferred(no_offset_view(P)) === P
    @test @inferred(no_offset_view(A)) === A
    a0 = reshape([1])
    @test @inferred(no_offset_view(a0)) === a0
    a0v = view(a0)
    @test @inferred(no_offset_view(a0v)) === a0v

    # generic fallback
    A = collect(reshape(1:12, 3, 4))
    N = NegativeArray(A)
    @test N[-3, -4] == 1
    V = no_offset_view(N)
    @test collect(V) == A
    A = reshape(view([5], 1, 1))
    @test no_offset_view(A) == A

    # bidirectional
    B = BidirectionalVector([1, 2, 3])
    pushfirst!(B, 0)
    OB = no_offset_view(B)
    @test no_offset_axes(OB, 1) == 1:4
    @test collect(OB) == 0:3

    # issue #198
    offax = axes(OffsetVector(1:10, -5), 1)
    noffax = no_offset_view(offax)
    @test noffax == -4:5
    @test axes(noffax, 1) == 1:10   # ideally covered by the above, but current it isn't
    @test isa(noffax, AbstractUnitRange)

    r = Base.OneTo(4)
    @test no_offset_view(r) isa typeof(r)

    # SubArrays
    A = reshape(1:12, 3, 4)
    V = view(A, IdentityUnitRange(2:3), IdentityUnitRange(2:3))
    if collect(V) == [5 8; 6 9]   # julia 1.0 has a bug here
        @test no_offset_view(V) == [5 8; 6 9]
    end
    V = view(A, IdentityUnitRange(2:3), 2)
    @test V != [5;6]
    if collect(V) == [5;6]
        @test no_offset_view(V) == [5;6]
    end
    O = OffsetArray(A, -1:1, 0:3)
    V = view(O, 0:1, 1:2)
    @test V == no_offset_view(V) == [5 8; 6 9]
    r1, r2 = IdOffsetRange(1:3, -2), IdentityUnitRange(2:3)
    V = view(O, r1, r2)
    @test V != collect(V)
    @test no_offset_view(V) == collect(V)
    V = @view O[:,:]
    @test IndexStyle(A) == IndexStyle(O) == IndexStyle(V) == IndexStyle(no_offset_view(V)) == IndexLinear()

    @testset "issue #375" begin
        arr = OffsetArray(reshape(1:15, 3, 5), 2, 3)
        arr_no_offset = no_offset_view(@view arr[:, 4])
        @test all(!Base.has_offset_axes, axes(arr_no_offset))
    end
end

@testset "no nesting" begin
    A = randn(2, 3)
    x = A[2, 2]
    O1 = OffsetArray(A, -1:0, -1:1)
    O2 = OffsetArray(O1, 0:1, 0:2)
    @test parent(O1) ≡ parent(O2)
    @test eltype(O1) ≡ eltype(O2)
    O2[1, 1] = x + 1            # just a sanity check
    @test A[2, 2] == x + 1
end

@testset "mutating functions for OffsetVector" begin
    # push!
    o = OffsetVector(Int[], -1)
    @test push!(o) === o
    @test no_offset_axes(o, 1) == 0:-1
    @test push!(o, 1) === o
    @test no_offset_axes(o, 1) == 0:0
    @test o[end] == 1
    @test push!(o, 2, 3) === o
    @test no_offset_axes(o, 1) == 0:2
    @test o[end-1:end] == [2, 3]
    # pop!
    o = OffsetVector([1, 2, 3], -1)
    @test pop!(o) == 3
    @test no_offset_axes(o, 1) == 0:1
    # append!
    o = OffsetVector([1, 2, 3], -1)
    append!(o, [4, 5])
    @test no_offset_axes(o, 1) == 0:4
    # empty!
    o = OffsetVector([1, 2, 3], -1)
    @test empty!(o) === o
    @test no_offset_axes(o, 1) == 0:-1
end

@testset "searchsorted (#85)" begin
    o = OffsetVector([1,3,4,5],-2)
    @test searchsortedfirst(o,-2) == -1
    @test searchsortedfirst(o, 1) == -1
    @test searchsortedfirst(o, 2) ==  0
    @test searchsortedfirst(o, 5) ==  2
    @test searchsortedfirst(o, 6) ==  3
    @test searchsortedlast(o, -2) == -2
    @test searchsortedlast(o,  1) == -1
    @test searchsortedlast(o,  2) == -1
    @test searchsortedlast(o,  5) ==  2
    @test searchsortedlast(o,  6) ==  2
    @test searchsorted(o, -2) == -1:-2
    @test searchsorted(o,  1) == -1:-1
    @test searchsorted(o,  2) ==  0:-1
    @test searchsorted(o,  5) ==  2:2
    @test searchsorted(o,  6) ==  3:2

    if VERSION > v"1.2"
        # OffsetVector of another offset vector
        v = OffsetVector(Base.IdentityUnitRange(4:10),-2)
        @test searchsortedfirst(v, first(v)-1) == firstindex(v)
        for i in axes(v,1)
            @test searchsortedfirst(v, v[i]) == i
        end
        @test searchsortedfirst(v, last(v)+1) == lastindex(v)+1
        @test searchsortedlast(v, first(v)-1) == firstindex(v)-1
        for i in axes(v,1)
            @test searchsortedlast(v, v[i]) == i
        end
        @test searchsortedlast(v, last(v)+1) == lastindex(v)
        @test searchsorted(v, first(v)-1) === firstindex(v) .+ (0:-1)
        for i in axes(v,1)
            @test searchsorted(v, v[i]) == i:i
        end
        @test searchsorted(v, last(v)+1) === lastindex(v) .+ (1:0)
    end

    v = OffsetVector{Float64, OffsetVector{Float64, Vector{Float64}}}(OffsetVector([2,2,3,3,3,4], 3), 4)
    @test searchsortedfirst(v, minimum(v)-1) == firstindex(v)
    for el in unique(v)
        @test searchsortedfirst(v, el) == findfirst(isequal(el), v)
    end
    @test searchsortedfirst(v, maximum(v)+1) == lastindex(v)+1

    @test searchsortedlast(v, minimum(v)-1) == firstindex(v)-1
    for el in unique(v)
        @test searchsortedlast(v, el) == findlast(isequal(el), v)
    end
    @test searchsortedlast(v, maximum(v)+1) == lastindex(v)

    @test searchsorted(v, minimum(v)-1) === firstindex(v) .+ (0:-1)
    for el in unique(v)
        @test searchsorted(v, el) == findfirst(isequal(el), v):findlast(isequal(el), v)
    end
    @test searchsorted(v, maximum(v)+1) === lastindex(v) .+ (1:0)

    soa = OffsetArray([2,2,3], typemax(Int)-3)
    @test searchsortedfirst(soa, 1) == firstindex(soa) == typemax(Int)-2
    @test searchsortedfirst(soa, 2) == firstindex(soa) == typemax(Int)-2
    @test searchsortedfirst(soa, 3) == lastindex(soa) == typemax(Int)

    soa = OffsetArray([2,2,3], typemin(Int))
    @test searchsortedlast(soa, 2) == firstindex(soa) + 1 == typemin(Int) + 2
    @test searchsortedlast(soa, 3) == lastindex(soa) == typemin(Int) + 3
    @test searchsortedlast(soa, 1) == typemin(Int)

    soa = OffsetArray([2,2,3], typemax(Int)-4)
    @test searchsorted(soa, 1) === firstindex(soa) .+ (0:-1)
    @test searchsorted(soa, 2) == firstindex(soa) .+ (0:1) == typemax(Int) .+ (-3:-2)
    @test searchsorted(soa, 3) == lastindex(soa) .+ (0:0) == typemax(Int) .+ (-1:-1)
    @test searchsorted(soa, 4) === lastindex(soa) .+ (1:0)

    soa = OffsetArray([2,2,3], typemax(Int)-3)
    @test searchsorted(soa, 1) === firstindex(soa) .+ (0:-1)
    @test searchsorted(soa, 2) == firstindex(soa) .+ (0:1) == typemax(Int) .+ (-2:-1)
    @test searchsorted(soa, 3) == lastindex(soa) .+ (0:0) == typemax(Int) .+ (0:0)
    @test searchsorted(soa, 4) === lastindex(soa) .+ (1:0)

    soa = OffsetArray([2,2,3], typemin(Int))
    @test searchsorted(soa, 1) === firstindex(soa) .+ (0:-1)
    @test searchsorted(soa, 2) == firstindex(soa) .+ (0:1) == typemin(Int) .+ (1:2)
    @test searchsorted(soa, 3) == lastindex(soa) .+ (0:0) == typemin(Int) .+ (3:3)
    @test searchsorted(soa, 4) === lastindex(soa) .+ (1:0)
end

@testset "Pointer" begin
    a = OffsetVector(collect(10:20), 9);
    @test 12 == a[12] == unsafe_load(pointer(a), 12 + (1 - firstindex(a))) == unsafe_load(pointer(a, 12))

    A = OffsetArray(reshape(collect(10:130), (11,11)), 9, 9);
    @test 21 == A[12] == unsafe_load(pointer(A), 12) == unsafe_load(pointer(A, 12))
    @test 61 == A[52] == unsafe_load(pointer(A), 52) == unsafe_load(pointer(A, 52))

    @test pointer(a) === pointer(parent(a))
    @test pointer(A) === pointer(parent(A))
    @test pointer(a, 12) === pointer(parent(a), 12 + (1 - firstindex(a)))
    @test pointer(A, 12) === pointer(parent(A), 12)
    @test pointer(a) === pointer(a, firstindex(a))
    @test pointer(A) === pointer(A, firstindex(A))
    if VERSION ≥ v"1.5"
        @test pointer(a') === pointer(parent(a))
        @test pointer(A') === pointer(parent(A))
        @test pointer(a', 5) === pointer(parent(a), 5)
        @test pointer(A', 15) === pointer(parent(A)', 15)
    end

    @test Base.cconvert(Ptr{eltype(A)}, A) == Base.cconvert(Ptr{eltype(A)}, parent(A))
end

struct Foo2
    o::OffsetArray{Float64,1,Array{Float64,1}}
end

@testset "convert" begin
    d = Diagonal([1,1,1])
    M = convert(Matrix{Float64}, d)
    od = OffsetArray(d, 1, 1)
    oM = convert(OffsetMatrix{Float64, Matrix{Float64}}, od)
    @test eltype(oM) == Float64
    @test typeof(parent(oM)) == Matrix{Float64}
    @test oM == od
    oM2 = convert(OffsetMatrix{Float64, Matrix{Float64}}, d)
    @test eltype(oM2) == Float64
    @test typeof(parent(oM2)) == Matrix{Float64}
    @test oM2 == d

    # issue 171
    O = OffsetArray(zeros(Int, 3), -1)
    F = Foo2(O)
    @test F.o == O

    a = [MMatrix{2,2}(1:4) for i = 1:2]
    oa = [OffsetArray(ai, 0, 0) for ai in a]
    b = ones(2,2)
    @test b * a == b * oa

    for a = [1:4, OffsetArray(ones(5))]
        for T in [OffsetArray, OffsetVector,
            OffsetArray{eltype(a)}, OffsetArray{Float32},
            OffsetVector{eltype(a)}, OffsetVector{Float32},
            OffsetVector{Float32, Vector{Float32}},
            OffsetVector{Float32, OffsetVector{Float32, Vector{Float32}}},
            OffsetVector{eltype(a), typeof(a)},
            ]

            @test convert(T, a) isa T
            @test convert(T, a) == a

            b = T(a)
            @test b isa T
            @test b == a

            b = T(a, 0)
            @test b isa T
            @test b == a

            b = T(a, axes(a))
            @test b isa T
            @test b == a
        end

        a2 = reshape(a, :, 1)
        for T in [OffsetArray{Float32}, OffsetMatrix{Float32}, OffsetArray{Float32, 2, Matrix{Float32}}]
            b = T(a2, 0, 0)
            @test b isa T
            @test b == a2

            b = T(a2, axes(a2))
            @test b isa T
            @test b == a2

            b = T(a2, 1, 1)
            @test no_offset_axes(b) == map((x,y) -> x .+ y, axes(a2), (1,1))

            b = T(a2)
            @test b isa T
            @test b == a2
        end
        a3 = reshape(a, :, 1, 1)
        for T in [OffsetArray{Float32}, OffsetArray{Float32, 3}, OffsetArray{Float32, 3, Array{Float32,3}}]
            b = T(a3, 0, 0, 0)
            @test b isa T
            @test b == a3

            b = T(a3, axes(a3))
            @test b isa T
            @test b == a3

            b = T(a3, 1, 1, 1)
            @test no_offset_axes(b) == map((x,y) -> x .+ y, axes(a3), (1,1,1))

            b = T(a3)
            @test b isa T
            @test b == a3
        end
    end

    a = OffsetArray(ones(2), 1)
    b = convert(OffsetArray, a)
    @test a === b
    b = convert(OffsetVector, a)
    @test a === b

    # test that non-Int offsets work correctly
    a = 1:4
    b1 = OffsetVector{Float64,Vector{Float64}}(a, 2)
    b2 = OffsetVector{Float64,Vector{Float64}}(a, big(2))
    @test b1 == b2

    a = OffsetArray(ones(2), 1)
    b1 = OffsetArray{Float64, 1, typeof(a)}(a, (-1,))
    b2 = OffsetArray{Float64, 1, typeof(a)}(a, (-big(1),))
    @test b1 == b2

    # test for custom offset arrays
    a = ZeroBasedRange(1:3)
    for T in [OffsetVector{Float64, UnitRange{Float64}}, OffsetVector{Int, Vector{Int}},
        OffsetVector{Float64,OffsetVector{Float64,UnitRange{Float64}}},
        OffsetArray{Int,1,OffsetArray{Int,1,UnitRange{Int}}},
        ]

        b = T(a)
        @test b isa T
        @test b == a

        b = T(a, 2:4)
        @test b isa T
        @test no_offset_axes(b, 1) == 2:4
        @test no_offset_view(b) == no_offset_view(a)

        b = T(a, 1)
        @test b isa T
        @test axes(b, 1) == 1:3
        @test no_offset_view(b) == no_offset_view(a)

        c = convert(T, a)
        @test c isa T
        @test c == a
    end

    # test using custom indices
    a = ones(2,2)
    for T in [OffsetMatrix{Int}, OffsetMatrix{Float64}, OffsetMatrix{Float64, Matrix{Float64}},
        OffsetMatrix{Int, Matrix{Int}}]

        b = T(a, ZeroBasedIndexing())
        @test b isa T
        @test no_offset_axes(b) == (0:1, 0:1)
    end

    # changing the number of dimensions is not permitted
    A = rand(2,2)
    @test_throws MethodError convert(OffsetArray{Float64, 3}, A)
    @test_throws MethodError convert(OffsetArray{Float64, 3, Array{Float64,3}}, A)
end

@testset "Conversion to AbstractArray{T}" begin
    r = 1:4
    T = Float64
    V = typeof(map(T, r))
    v = OffsetVector(r)
    @test OffsetArray{T}(v) isa OffsetVector{T,V}
    @test AbstractArray{T}(v) isa OffsetVector{T,V}
    @test AbstractVector{T}(v) isa OffsetVector{T,V}
    @test convert(AbstractVector{T}, v) isa OffsetVector{T,V}
    @test convert(AbstractArray{T}, v) isa OffsetVector{T,V}
    @test axes(OffsetArray{T}(v)) === axes(v)
    @test axes(AbstractArray{T}(v)) === axes(v)
    @test axes(AbstractVector{T}(v)) === axes(v)
    @test axes(convert(AbstractVector{T}, v)) === axes(v)
    @test axes(convert(AbstractArray{T}, v)) == axes(v)

    A = SMatrix{2,2}(1, 0, 0, 1)
    TA = typeof(map(T, A))
    OA = OffsetMatrix(A, 3:4, 5:6)
    @test OffsetArray{T}(OA) isa OffsetMatrix{T,TA}
    @test AbstractArray{T}(OA) isa OffsetMatrix{T,TA}
    @test AbstractMatrix{T}(OA) isa OffsetMatrix{T,TA}
    @test convert(AbstractMatrix{T}, OA) isa OffsetMatrix{T,TA}
    @test convert(AbstractArray{T}, OA) isa OffsetMatrix{T,TA}
    @test axes(OffsetArray{T}(OA)) === axes(OA)
    @test axes(AbstractArray{T}(OA)) === axes(OA)
    @test axes(AbstractMatrix{T}(OA)) === axes(OA)
    @test axes(convert(AbstractMatrix{T}, OA)) === axes(OA)
    @test axes(convert(AbstractArray{T}, OA)) === axes(OA)
end

@testset "unsafe_wrap" begin
    p = Ptr{UInt16}(Libc.malloc(2*3*4*2))
    @test unsafe_wrap(OffsetArray, p, 2, 3, 4) isa OffsetArray{UInt16, 3}
    @test unsafe_wrap(OffsetArray, p, (2, 3, 4)) isa OffsetArray{UInt16, 3}
    @test unsafe_wrap(OffsetVector, p, 2*3*4) isa OffsetVector{UInt16}
    @test unsafe_wrap(OffsetMatrix, p, 2*3, 4) isa OffsetMatrix{UInt16}
    @test unsafe_wrap(OffsetArray{UInt16}, p, 2, 3, 4) isa OffsetArray{UInt16, 3}
    @test unsafe_wrap(OffsetArray{UInt16}, p, (2, 3, 4)) isa OffsetArray{UInt16, 3}
    @test unsafe_wrap(OffsetVector{UInt16}, p, 2*3*4) isa OffsetVector{UInt16}
    @test unsafe_wrap(OffsetMatrix{UInt16}, p, 2*3, 4) isa OffsetMatrix{UInt16}
    p = Ptr{UInt8}(p)
    @test unsafe_wrap(OffsetArray, p, 2:3, 3:5, 4:7) isa OffsetArray{UInt8, 3}
    @test unsafe_wrap(OffsetArray, p, (2:3, 3:5, 4:7)) isa OffsetArray{UInt8, 3}
    @test unsafe_wrap(OffsetVector, p, 1:(2*3*4) .- 1) isa OffsetVector{UInt8}
    @test unsafe_wrap(OffsetMatrix, p, 1:(2*3) .+ 6, 4:7) isa OffsetMatrix{UInt8}
    @test unsafe_wrap(OffsetMatrix, p, -5:5, Base.OneTo(3); own = true) isa OffsetMatrix{UInt8}
end

end