using OffsetArrays
using OffsetArrays: IdentityUnitRange, no_offset_view
using OffsetArrays: IdOffsetRange
using Test, Aqua, Documenter
using LinearAlgebra
using DelimitedFiles
using CatIndices: BidirectionalVector
using EllipsisNotation
using Adapt
using StaticArrays

DocMeta.setdocmeta!(OffsetArrays, :DocTestSetup, :(using OffsetArrays); recursive=true)

# https://github.com/JuliaLang/julia/pull/29440
if VERSION < v"1.1.0-DEV.389"
    Base.:(:)(I::CartesianIndex{N}, J::CartesianIndex{N}) where N =
        CartesianIndices(map((i,j) -> i:j, Tuple(I), Tuple(J)))
end

# Custom index types
struct ZeroBasedIndexing end
struct NewColon end
struct TupleOfRanges{N}
    x ::NTuple{N, UnitRange{Int}}
end

# Useful for testing indexing
struct ZeroBasedRange{T,A<:AbstractRange{T}} <: AbstractRange{T}
    a :: A
    function ZeroBasedRange(a::AbstractRange{T}) where {T}
        @assert !Base.has_offset_axes(a)
        new{T, typeof(a)}(a)
    end
end

struct ZeroBasedUnitRange{T,A<:AbstractUnitRange{T}} <: AbstractUnitRange{T}
    a :: A
    function ZeroBasedUnitRange(a::AbstractUnitRange{T}) where {T}
        @assert !Base.has_offset_axes(a)
        new{T, typeof(a)}(a)
    end
end

for Z in [:ZeroBasedRange, :ZeroBasedUnitRange]
    @eval Base.parent(A::$Z) = A.a
    @eval Base.first(A::$Z) = first(A.a)
    @eval Base.length(A::$Z) = length(A.a)
    @eval Base.last(A::$Z) = last(A.a)
    @eval Base.size(A::$Z) = size(A.a)
    @eval Base.axes(A::$Z) = map(x -> IdentityUnitRange(0:x-1), size(A.a))
    @eval Base.getindex(A::$Z, i::Int) = A.a[i + 1]
    @eval Base.step(A::$Z) = step(A.a)
    @eval OffsetArrays.no_offset_view(A::$Z) = A.a
    @eval function Base.show(io::IO, A::$Z)
        show(io, A.a)
        print(io, " with indices $(axes(A,1))")
    end

    for R in [:AbstractRange, :AbstractUnitRange, :StepRange]
        @eval @inline function Base.getindex(A::$Z, r::$R{<:Integer})
            @boundscheck checkbounds(A, r)
            OffsetArray(A.a[r .+ 1], axes(r))
        end
    end
    for R in [:UnitRange, :StepRange, :StepRangeLen, :LinRange]
        @eval @inline function Base.getindex(A::$R, r::$Z)
            @boundscheck checkbounds(A, r)
            OffsetArray(A[r.a], axes(r))
        end
    end
    @eval @inline function Base.getindex(A::StepRangeLen{<:Any,<:Base.TwicePrecision,<:Base.TwicePrecision}, r::$Z)
        @boundscheck checkbounds(A, r)
        OffsetArray(A[r.a], axes(r))
    end
end

function same_value(r1, r2)
    length(r1) == length(r2) || return false
    for (v1, v2) in zip(r1, r2)
        v1 == v2 || return false
    end
    return true
end

@testset "Project meta quality checks" begin
    # Not checking compat section for test-only dependencies
    Aqua.test_all(OffsetArrays; project_extras=true, deps_compat=true, stale_deps=true, project_toml_formatting=true)
    if VERSION >= v"1.2"
        doctest(OffsetArrays, manual = false)
    end
end

@testset "IdOffsetRange" begin

    function check_indexed_by(r, rindx)
        for i in rindx
            r[i]
        end
        @test_throws BoundsError r[minimum(rindx)-1]
        @test_throws BoundsError r[maximum(rindx)+1]
        return nothing
    end

    ro = OffsetArrays.IdOffsetRange(Base.OneTo(3))
    rs = OffsetArrays.IdOffsetRange(3:5, -2)
    @test typeof(ro) !== typeof(rs)
    @test same_value(ro, 1:3)
    check_indexed_by(ro, 1:3)
    @test same_value(rs, 1:3)
    check_indexed_by(rs, -1:1)
    @test firstindex(ro) == 1
    @test lastindex(ro) == 3
    @test firstindex(rs) == -1
    @test lastindex(rs) == 1
    @test @inferred(typeof(ro)(ro)) === ro
    @test @inferred(OffsetArrays.IdOffsetRange{Int}(ro))   === ro
    @test @inferred(OffsetArrays.IdOffsetRange{Int16}(ro)) === OffsetArrays.IdOffsetRange(Base.OneTo(Int16(3)))
    @test @inferred(OffsetArrays.IdOffsetRange(ro))        === ro
    @test parent(ro) === ro.parent
    @test parent(rs) === rs.parent
    # construction/coercion preserves the values, altering the axes if needed
    r2 = @inferred(typeof(rs)(ro))
    @test typeof(r2) === typeof(rs)
    @test same_value(ro, 1:3)
    check_indexed_by(ro, 1:3)
    r2 = @inferred(typeof(ro)(rs))
    @test typeof(r2) === typeof(ro)
    @test same_value(r2, 1:3)
    check_indexed_by(r2, 1:3)
    # check the example in the comments
    r = OffsetArrays.IdOffsetRange{Int,UnitRange{Int}}(3:4)
    @test same_value(r, 3:4)
    check_indexed_by(r, 1:2)
    r = OffsetArrays.IdOffsetRange{Int,Base.OneTo{Int}}(3:4)
    @test same_value(r, 3:4)
    check_indexed_by(r, 3:4)
    r = OffsetArrays.IdOffsetRange{Int,Base.OneTo{Int}}(3:4, -2)
    @test same_value(r, 1:2)
    check_indexed_by(r, 1:2)

    r = OffsetArrays.IdOffsetRange{Int32, Base.OneTo{Int32}}(Base.OneTo(Int64(2)), 3)
    @test same_value(r, 4:5)
    check_indexed_by(r, 4:5)

    r = IdOffsetRange{Int, UnitRange{Int}}(IdOffsetRange(3:5, 2), 2)
    @test typeof(r) == IdOffsetRange{Int, UnitRange{Int}}
    @test same_value(r, 7:9)
    check_indexed_by(r, 5:7)

    r = IdOffsetRange{Int, Base.OneTo{Int}}(IdOffsetRange(Base.OneTo(3), 1), 1)
    @test typeof(r) == IdOffsetRange{Int,Base.OneTo{Int}}
    @test same_value(r, 3:5)
    check_indexed_by(r, 3:5)

    rp = Base.OneTo(3)
    r = IdOffsetRange(rp)
    r2 = IdOffsetRange{Int,typeof(r)}(r, 1)
    @test same_value(r2, 2:4)
    check_indexed_by(r2, 2:4)

    r2 = IdOffsetRange{Int32,IdOffsetRange{Int32,Base.OneTo{Int32}}}(r, 1)
    @test typeof(r2) == IdOffsetRange{Int32,IdOffsetRange{Int32,Base.OneTo{Int32}}}
    @test same_value(r2, 2:4)
    check_indexed_by(r2, 2:4)

    # Constructor that's round-trippable with `show`
    rrt = IdOffsetRange(values=7:9, indices=-1:1)
    @test same_value(rrt, 7:9)
    check_indexed_by(rrt, -1:1)
    @test_throws ArgumentError IdOffsetRange(values=7:9, indices=-1:2)
    @test_throws ArgumentError IdOffsetRange(values=7:9, indices=-1:0)
    @test_throws TypeError IdOffsetRange(values=7:9, indices=-1)
    @test_throws UndefKeywordError IdOffsetRange(values=7:9)
    @test_throws UndefKeywordError IdOffsetRange(indices=-1:1)
    @test_throws MethodError IdOffsetRange(7:9, indices=-1:1)
    @test_throws MethodError IdOffsetRange(-1:1, values=7:9)

    # conversion preserves both the values and the axes, throwing an error if this is not possible
    @test @inferred(oftype(ro, ro)) === ro
    @test @inferred(convert(OffsetArrays.IdOffsetRange{Int}, ro)) === ro
    @test @inferred(convert(OffsetArrays.IdOffsetRange{Int}, rs)) === rs
    @test @inferred(convert(OffsetArrays.IdOffsetRange{Int16}, ro)) === OffsetArrays.IdOffsetRange(Base.OneTo(Int16(3)))
    r2 = @inferred(oftype(rs, ro))
    @test typeof(r2) === typeof(rs)
    @test same_value(r2, 1:3)
    check_indexed_by(r2, 1:3)
    # These two broken tests can be fixed by uncommenting the `convert` definitions
    # in axes.jl, but unfortunately Julia may not quite be ready for this. (E.g. `reinterpretarray.jl`)
    @test_broken try oftype(ro, rs); false catch err true end  # replace with line below
    # @test_throws ArgumentError oftype(ro, rs)
    @test @inferred(oftype(ro, Base.OneTo(2))) === OffsetArrays.IdOffsetRange(Base.OneTo(2))
    @test @inferred(oftype(ro, 1:2)) === OffsetArrays.IdOffsetRange(Base.OneTo(2))
    @test_broken try oftype(ro, 3:4); false catch err true end
    # @test_throws ArgumentError oftype(ro, 3:4)

    # broadcasting behavior with scalars (issue #104)
    r3 = (1 .+ OffsetArrays.IdOffsetRange(3:5, -1) .+ 1) .- 1
    @test same_value(r3, 3:5)
    check_indexed_by(r3, 0:2)

    @testset "Idempotent indexing" begin
        @testset "Indexing into an IdOffsetRange" begin
            r = OffsetArrays.IdOffsetRange(3:5, -1)
            # Indexing with IdentityUnitRange
            s = IdentityUnitRange(0:2)
            @test axes(r[s]) == axes(s)
            for i in eachindex(s)
                @test r[s[i]] == r[s][i]
            end

            # Indexing with IdOffsetRange
            s = OffsetArrays.IdOffsetRange(-4:-2, 4)
            @test axes(r[s]) == axes(s)
            for i in eachindex(s)
                @test r[s[i]] == r[s][i]
            end

            # Indexing with UnitRange
            s = 0:2
            @test axes(r[s]) == axes(s)
            for i in eachindex(s)
                @test r[s[i]] == r[s][i]
            end
        end
        @testset "Indexing using an IdOffsetRange" begin
            r = OffsetArrays.IdOffsetRange(3:5, -1)
            # Indexing into an IdentityUnitRange
            s = IdentityUnitRange(-1:5)
            @test axes(s[r]) == axes(r)
            for i in eachindex(r)
                @test s[r[i]] == s[r][i]
            end

            # Indexing into an UnitRange
            s = -3:6
            @test axes(s[r]) == axes(r)
            for i in eachindex(r)
                @test s[r[i]] == s[r][i]
            end
        end
    end

    # Test reduced index
    rred = Base.reduced_index(r)
    @test typeof(rred) == typeof(r)
    @test length(rred) == 1
    @test first(rred) == first(r)
end

# used in testing the constructor
struct WeirdInteger{T} <: Integer
    x :: T
end
# assume that it doesn't behave as expected
Base.convert(::Type{Int}, a::WeirdInteger) = a

@testset "Constructors" begin
    @testset "Single-entry arrays in dims 0:5" begin
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
        @test_nowarn OffsetArray{Float64, 1, typeof(ao)}(ao, (-1, ))
        @test_throws OverflowError OffsetArray{Float64, 1, typeof(ao)}(ao, (-2, )) # inner Constructor
        @test_throws OverflowError OffsetArray(ao, (-2, )) # convinient constructor accumulate offsets

        @testset "OffsetRange" begin
            local r = 1:100
            local a = OffsetVector(r, 4)
            @test first(r) in a
            @test !(last(r) + 1 in a)
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
        @test axes(OffsetArray(a, :, CartesianIndices((3:4, 2:3)))) == (1:2, 3:4, 2:3)
        @test axes(OffsetArray(a, 10:11, CartesianIndices((3:4, 2:3)) )) == (10:11, 3:4, 2:3)
        @test axes(OffsetArray(a, CartesianIndices((3:4, 2:3)), :)) == (3:4, 2:3, 1:2)
        @test axes(OffsetArray(a, CartesianIndices((3:4, 2:3)), 10:11)) == (3:4, 2:3, 10:11)
        @test axes(OffsetArray(a, :, :, CartesianIndices((3:4,)) )) == (1:2, 1:2, 3:4)
        @test axes(OffsetArray(a, 10:11, :, CartesianIndices((3:4,)) )) == (10:11, 1:2, 3:4)
        @test axes(OffsetArray(a, 10:11, 2:3, CartesianIndices((3:4,)) )) == (10:11, 2:3, 3:4)

        # ignore empty CartesianIndices
        @test OffsetArray(a, CartesianIndices(()), 0:1, :, 2:3) == OffsetArray(a, 0:1, :, 2:3)
        @test OffsetArray(a, 0:1, CartesianIndices(()), :, 2:3) == OffsetArray(a, 0:1, :, 2:3)
        @test OffsetArray(a, 0:1, :,  CartesianIndices(()), 2:3) == OffsetArray(a, 0:1, :, 2:3)
        @test OffsetArray(a, 0:1, :, 2:3, CartesianIndices(())) == OffsetArray(a, 0:1, :, 2:3)

        indices = (-1:1, -7:7, -1:2, -5:5, -1:1, -3:3, -2:2, -1:1)
        y = OffsetArray{Float64}(undef, indices...);
        @test axes(y) === axes(OffsetArray{Float64}(undef, indices))
        @test axes(y) === axes(OffsetArray{Float64, length(indices)}(undef, indices...))
        @test axes(y) == (-1:1, -7:7, -1:2, -5:5, -1:1, -3:3, -2:2, -1:1)
        @test eltype(y) === Float64

        @test_throws ArgumentError OffsetArray{Float64, 2}(undef, indices)
        @test_throws ArgumentError OffsetArray(y, indices[1:2])

        @test ndims(OffsetArray(zeros(), ())) == 0
        @test Base.axes1(OffsetArray(zeros(), ())) === OffsetArrays.IdOffsetRange(Base.OneTo(1))

        @testset "convenience constructors" begin
            ax = (2:3, 4:5)

            for f in (zeros, ones)
                a = f(Float64, ax)
                @test axes(a) == ax
                @test eltype(a) == Float64
            end

            for f in (trues, falses)
                a = f(ax)
                @test axes(a) == ax
                @test eltype(a) == Bool
            end
        end

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

                @test axes(OffsetArray(v, .., 2:6)) == (2:6, )
                @test OffsetArray(v, .., 2:6) == OffsetArray(v, 2:6)
                @test axes(OffsetVector(v, .., 2:6)) == (2:6, )
                @test OffsetVector(v, .., 2:6) == OffsetVector(v, 2:6)
            end
            @testset "Matrix" begin
                m = rand(2, 2)
                @test axes(OffsetArray(m, ..)) == axes(m)
                @test OffsetArray(m, ..) == OffsetArray(m, :, :)
                @test axes(OffsetMatrix(m, ..)) == axes(m)
                @test OffsetMatrix(m, ..) == OffsetMatrix(m, :, :)

                @test axes(OffsetArray(m, .., 2:3)) == (axes(m, 1), 2:3)
                @test OffsetArray(m, .., 2:3) == OffsetArray(m, :, 2:3)
                @test axes(OffsetMatrix(m, .., 2:3)) == (axes(m, 1), 2:3)
                @test OffsetMatrix(m, .., 2:3) == OffsetMatrix(m, :, 2:3)

                @test axes(OffsetArray(m, .., 2:3, 3:4)) == (2:3, 3:4)
                @test OffsetArray(m, .., 2:3, 3:4) == OffsetArray(m, 2:3, 3:4)
                @test axes(OffsetMatrix(m, .., 2:3, 3:4)) == (2:3, 3:4)
                @test OffsetMatrix(m, .., 2:3, 3:4) == OffsetMatrix(m, 2:3, 3:4)
            end
            @testset "3D Array" begin
                a = rand(2, 2, 2)
                @test axes(OffsetArray(a, ..)) == axes(a)
                @test OffsetArray(a, ..) == OffsetArray(a, :, :, :)

                @test axes(OffsetArray(a, .., 2:3)) == (axes(a)[1:2]..., 2:3)
                @test OffsetArray(a, .., 2:3) == OffsetArray(a, :, :, 2:3)

                @test axes(OffsetArray(a, .., 2:3, 3:4)) == (axes(a, 1), 2:3, 3:4)
                @test OffsetArray(a, .., 2:3, 3:4) == OffsetArray(a, :, 2:3, 3:4)

                @test axes(OffsetArray(a, 2:3, .., 3:4)) == (2:3, axes(a, 2), 3:4)
                @test OffsetArray(a, 2:3, .., 3:4) == OffsetArray(a, 2:3, :, 3:4)

                @test axes(OffsetArray(a, .., 4:5, 2:3, 3:4)) == (4:5, 2:3, 3:4)
                @test OffsetArray(a, .., 4:5, 2:3, 3:4) == OffsetArray(a, 4:5, 2:3, 3:4)
            end
        end
        @testset "ZeroBasedIndexing" begin
            Base.to_indices(A, inds, ::Tuple{ZeroBasedIndexing}) = map(x -> 0:length(x) - 1, inds)

            a = zeros(3,3)
            oa = OffsetArray(a, ZeroBasedIndexing())
            @test axes(oa) == (0:2, 0:2)
        end
        @testset "TupleOfRanges" begin
            Base.to_indices(A, inds, t::Tuple{TupleOfRanges{N}}) where {N} = t
            OffsetArrays.AxisConversionStyle(::Type{TupleOfRanges{N}}) where {N} =
                OffsetArrays.TupleOfRanges()

            Base.convert(::Type{Tuple{Vararg{AbstractUnitRange{Int}}}}, t::TupleOfRanges) = t.x

            a = zeros(3,3)
            inds = TupleOfRanges((3:5, 2:4))
            oa = OffsetArray(a, inds)
            @test axes(oa) == inds.x
        end
        @testset "NewColon" begin
            Base.to_indices(A, inds, t::Tuple{NewColon,Vararg{Any}}) =
                (_uncolon(inds, t), to_indices(A, Base.tail(inds), Base.tail(t))...)

            _uncolon(inds::Tuple{}, I::Tuple{NewColon, Vararg{Any}}) = OneTo(1)
            _uncolon(inds::Tuple, I::Tuple{NewColon, Vararg{Any}}) = inds[1]

            a = zeros(3, 3)
            oa = OffsetArray(a, (NewColon(), 2:4))
            @test axes(oa) == (axes(a,1), 2:4)
        end
    end

    @testset "Offset range construction" begin
        r = -2:5
        for AT in Any[OffsetArray, OffsetVector]
            y = AT(r, r)
            @test axes(y) == (r,)
            @test step(y) == step(r)
            y = AT(r, (r,))
            @test axes(y) == (r,)
            y = AT(r, CartesianIndices((r, )))
            @test axes(y) == (r, )
        end
    end
end

@testset "Axes supplied to constructor correspond to final result" begin
    # Ref https://github.com/JuliaArrays/OffsetArrays.jl/pull/65#issuecomment-457181268
    B = BidirectionalVector([1, 2, 3], -2)
    A = OffsetArray(B, -1:1)
    @test axes(A) == (IdentityUnitRange(-1:1),)
end

@testset "Traits" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))                   # IndexLinear
    S = OffsetArray(view(A0, 1:2, 1:2), (-1,2))   # IndexCartesian
    @test axes(A) === axes(S)
    @test axes(A) == axes(S) == (0:1, 3:4)
    @test axes(A, 1) === OffsetArrays.IdOffsetRange(Base.OneTo(2), -1)
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
    A = ones(5:6)
    @test eachindex(IndexLinear(), A) === axes(A, 1)
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

_comp(::Type{<:Integer}) = ==
_comp(::Type{<:Real}) = ≈
function test_indexing_axes_and_vals(r1, r2)
    r12 = r1[r2]
    op = _comp(eltype(r1))
    if axes(r12, 1) != axes(r2, 1)
        @show r1 r2 r12 axes(r12, 1) axes(r2, 1)
    end
    @test axes(r12, 1) == axes(r2, 1)
    if axes(r12, 1) == axes(r2, 1)
        @test op(first(r12), r1[first(r2)])
        @test op(last(r12), r1[last(r2)])
        for i in eachindex(r2)
            @test op(r12[i], r1[r2[i]])
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

    r1 = OffsetArray(IdentityUnitRange(100:1000), 3)
    r2 = r1[:]
    @test r2 == r1

    for r1 in Any[
        # AbstractArrays
        OffsetArray(10:1000, 0), # 1-based index
        OffsetArray(10:3:1000, 3), # offset index
        OffsetArray(10.0:3:1000.0, 0), # 1-based index
        OffsetArray(10.0:3:1000.0, 3), # offset index
        OffsetArray(IdOffsetRange(10:1000, 1), -1), # 1-based index
        OffsetArray(IdOffsetRange(10:1000, 1), 3), # offset index
        OffsetArray(IdOffsetRange(IdOffsetRange(10:1000, -4), 1), 3), # 1-based index
        OffsetArray(IdOffsetRange(IdOffsetRange(10:1000, -1), 1), 3), # offset index

        # AbstractRanges
        1:1000,
        UnitRange(1.0, 1000.0),
        1:3:1000,
        1.0:3.0:1000.0,
        StepRangeLen(Float64(1), Float64(1000), 1000),
        LinRange(1, 1000, 1000),
        IdOffsetRange(ZeroBasedUnitRange(1:1000), 1), # 1-based index
        IdOffsetRange(ZeroBasedUnitRange(1:1000), 2), # offset index
        ZeroBasedUnitRange(1:1000), # offset range
        ZeroBasedRange(1:1000), # offset range
        ZeroBasedRange(1:1:1000), # offset range
        ]

        # AbstractArrays with 1-based indices
        for r2 in Any[
            OffsetArray(5:80, 0),
            OffsetArray(5:2:80, 0),
            OffsetArray(IdentityUnitRange(5:80), -4),
            OffsetArray(IdOffsetRange(5:80), 0),
            ]

            test_indexing_axes_and_vals(r1, r2)
        end

        # AbstractRanges with 1-based indices
        for r2 in Any[
            5:80,
            5:2:80,
            IdOffsetRange(5:80),
            IdOffsetRange(ZeroBasedUnitRange(4:79), 1),
            ]

            test_indexing_axes_and_vals(r1, r2)
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
        1:1:99,
        1.0:1.0:99.0,
        StepRangeLen(Float64(1), Float64(99), 99),
        LinRange(1, 99, 99),
        ]

        for r2 in Any[
            IdentityUnitRange(Base.OneTo(10)),
            Base.Slice(Base.OneTo(10)),
            IdOffsetRange(Base.OneTo(10)),
            ]

            test_indexing_axes_and_vals(r1, r2)
        end
    end
end

@testset "Vector indexing with offset ranges" begin
    r = OffsetArray(8:10, -1:1)
    r1 = r[0:1]
    @test r1 === 9:10
    r1 = (8:10)[OffsetArray(1:2, -5:-4)]
    @test axes(r1) == (IdentityUnitRange(-5:-4),)
    @test parent(r1) === 8:9
    r1 = OffsetArray(8:10, -1:1)[OffsetArray(0:1, -5:-4)]
    @test axes(r1) == (IdentityUnitRange(-5:-4),)
    @test parent(r1) === 9:10

    a = OffsetVector(3:4, 10:11)
    ax = OffsetArrays.IdOffsetRange(5:6, 5)
    @test axes(a[ax]) == axes(ax)
    for i in axes(ax,1)
        @test a[ax[i]] == a[ax][i]
    end

    ax = IdentityUnitRange(10:11)
    @test axes(a[ax]) == axes(ax)
    for i in axes(ax,1)
        @test a[ax[i]] == a[ax][i]
    end

    for r1 in Any[
        # AbstractArrays
        OffsetArray(10:1000, 0), # 1-based index
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
        1:1000,
        UnitRange(1.0, 1000.0),
        1:2:2000,
        1.0:2.0:2000.0,
        StepRangeLen(Float64(1), Float64(1000), 1000),
        LinRange(1.0, 2000.0, 2000),
        IdOffsetRange(1:1000, 0), # 1-based index
        IdOffsetRange(ZeroBasedUnitRange(1:1000), 1), # 1-based index
        IdOffsetRange(ZeroBasedUnitRange(1:1000), 2), # offset index
        IdentityUnitRange(ZeroBasedUnitRange(1:1000)), # 1-based index
        ZeroBasedUnitRange(1:1000), # offset index
        ZeroBasedRange(1:1000), # offset index
        ZeroBasedRange(1:1:1000), # offset index
        ZeroBasedUnitRange(IdentityUnitRange(1:1000)), # offset index
        ]

        # AbstractArrays with offset axes
        for r2 in Any[OffsetArray(5:80, 40),
            OffsetArray(5:2:80, 40),
            OffsetArray(IdentityUnitRange(5:80), 2),
            OffsetArray(IdOffsetRange(5:80, 1), 3),
            OffsetArray(IdOffsetRange(IdOffsetRange(5:80, 4), 1), 3),
            OffsetArray(IdOffsetRange(IdentityUnitRange(5:80), 1), 3),
            OffsetArray(IdentityUnitRange(IdOffsetRange(5:80, 1)), 3),
            ]

            test_indexing_axes_and_vals(r1, r2)
        end

        # AbstractRanges with offset axes
        for r2 in Any[IdOffsetRange(5:80, 1),
            IdentityUnitRange(5:80),
            IdOffsetRange(Base.OneTo(9), 4),
            IdOffsetRange(IdOffsetRange(5:80, 2), 1),
            IdOffsetRange(IdOffsetRange(IdOffsetRange(5:80, -1), 2), 1),
            IdentityUnitRange(IdOffsetRange(1:10, 5)),
            IdOffsetRange(IdentityUnitRange(15:20), -2),
            ZeroBasedUnitRange(5:80),
            ZeroBasedRange(5:80),
            ZeroBasedRange(5:2:80),
            ]

            test_indexing_axes_and_vals(r1, r2)
        end
    end
end

@testset "LinearIndexing" begin
    r = OffsetArray(ZeroBasedRange(3:4), 1);
    @test LinearIndices(r) == axes(r,1)
    r = OffsetArray(ZeroBasedRange(3:4), 2);
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
    @test axes(y) == (r,)
    @test step(y) == step(s)

    a = OffsetVector(3:4, 10:11)
    ax = OffsetArrays.IdOffsetRange(5:6, 5)
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
    @test axes(S) == OffsetArrays.IdOffsetRange.((0:1, 3:4))
    # issue 100
    S = view(A, axes(A, 1), 3)
    @test S == A[:, 3]
    @test S[0] == 1
    @test S[1] == 2
    @test_throws BoundsError S[length(S)]
    @test axes(S) == (OffsetArrays.IdOffsetRange(0:1), )
    # issue 100
    S = view(A, 1, axes(A, 2))
    @test S == A[1, :]
    @test S[3] == 2
    @test S[4] == 4
    @test_throws BoundsError S[1]
    @test axes(S) == (OffsetArrays.IdOffsetRange(3:4), )

    # issue 133
    r = OffsetArrays.IdOffsetRange(1:2, -1)
    v1 = view(A, r, 3)
    @test v1[0] == 1
    @test v1[1] == 2
    @test axes(v1, 1) == axes(r, 1)
    v2 = view(A, UnitRange(r), 3)
    for (indflat, indoffset) in enumerate(r)
        @test v1[indoffset] == v2[indflat]
    end

    # issue 133
    r = OffsetArrays.IdOffsetRange(1:2, 2)
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
    r = OffsetArrays.IdOffsetRange(Base.OneTo(3), 5)
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
    @test axes(S) == (OffsetArrays.IdOffsetRange(0:1), Base.OneTo(2), OffsetArrays.IdOffsetRange(2:5))

    # issue #186
    a = reshape(1:12, 3, 4)
    r = OffsetArrays.IdOffsetRange(3:4)
    av = view(a, :, r)
    @test av == a[:, 3:4]
    @test axes(av) == (axes(a,1), axes(r,1))
    r = OffsetArrays.IdOffsetRange(1:2,2)
    av = view(a, :, r)
    @test no_offset_view(av) == a[:, 3:4]
    @test axes(av) == (axes(a,1), axes(r,1))
    r = OffsetArrays.IdOffsetRange(2:3)
    av1d = view(a, r, 3)
    @test av1d == a[2:3, 3]
    @test axes(av1d) == (axes(r,1),)
    r = OffsetArrays.IdOffsetRange(Base.OneTo(2), 1)
    av1d = view(a, r, 3)
    @test no_offset_view(av1d) == a[2:3, 3]
    @test axes(av1d) == (axes(r,1),)

    # fix IdOffsetRange(::IdOffsetRange, offset) nesting from #178
    b = 1:20
    bov = OffsetArray(view(b, 3:4), 3:4)
    c = @view b[bov]
    @test same_value(c, 3:4)
    @test axes(c,1) == 3:4
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
            av_nooffset = OffsetArrays.no_offset_view(av)
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

    a = OffsetArray([1 2; 3 4], -1:0, 5:6)
    shownsz = VERSION >= v"1.2.0-DEV.229" ? Base.dims2string(size(a))*' ' : ""
    @test summary(a) == "$(shownsz)OffsetArray(::$(typeof(parent(a))), -1:0, 5:6) with eltype $(Int) with indices -1:0×5:6"
    shownsz = VERSION >= v"1.2.0-DEV.229" ? Base.dims2string(size(view(a, :, 5)))*' ' : ""
    @test summary(view(a, :, 5)) == "$(shownsz)view(OffsetArray(::$(typeof(parent(a))), -1:0, 5:6), :, 5) with eltype $(Int) with indices -1:0"
    a = OffsetArray(reshape([1]))
    @test summary(a) == "0-dimensional OffsetArray(::$(typeof(parent(a)))) with eltype $(Int)"

    a = OffsetArray([1 2; 3 4], -1:0, 5:6)
    io = IOBuffer()
    show(io, axes(a, 1))
    @test String(take!(io)) == "IdOffsetRange(values=-1:0, indices=-1:0)"   # not qualified because of the using OffsetArrays: IdOffsetRange at top
    show(io, axes(a, 2))
    @test String(take!(io)) == "IdOffsetRange(values=5:6, indices=5:6)"
    rrtable = IdOffsetRange(values=7:9, indices=-1:1)
    rrted = eval(Meta.parse(string(rrtable)))
    @test pairs(rrtable) == pairs(rrted)

    @test Base.inds2string(axes(a)) == Base.inds2string(map(UnitRange, axes(a)))

    show(io, OffsetArray(3:5, 0:2))
    @test String(take!(io)) == "3:5 with indices 0:2"

    show(io, MIME"text/plain"(), OffsetArray(3:5, 0:2))
    @test String(take!(io)) == "3:5 with indices 0:2"

    # issue #198
    for r in Any[axes(OffsetVector(1:10, -5), 1), 1:1:2, 1.0:1.0:2.0, 1:-1:-5]
        a = OffsetVector(r, 5)
        show(io, a)
        @test String(take!(io)) == "$r with indices $(UnitRange(axes(a,1)))"
    end

    d = Diagonal([1,2,3])
    Base.print_array(io, d)
    s1 = String(take!(io))
    od = OffsetArray(d, -1:1, 3:5)
    Base.print_array(io, od)
    s2 = String(take!(io))
    @test s1 == s2

    @test Base.replace_in_print_matrix(od, -1, 3, " ") == Base.replace_in_print_matrix(d, 1, 1, " ")
    @test Base.replace_in_print_matrix(od, -1, 4, " ") == Base.replace_in_print_matrix(d, 1, 2, " ")

    v = rand(3)
    ov = OffsetArray(v, (-2,))
    @test Base.replace_in_print_matrix(ov, -1, 1, " ") == Base.replace_in_print_matrix(v, 1, 1, " ")

    # Avoid returning the value of toplevel if it is false
    # showarg should only print values, it shouldn't return anything
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
    B = similar(parent(A), (-3:3,1:4))
    @test isa(B, OffsetArray{Int,2})
    @test axes(B) == IdentityUnitRange.((-3:3, 1:4))
    @test isa([x for x in [1,2,3]], Vector{Int})
    @test similar(Array{Int}, (0:0, 0:0)) isa OffsetArray{Int, 2}
    @test similar(Array{Int}, (1, 1)) isa Matrix{Int}
    @test similar(Array{Int}, (Base.OneTo(1), Base.OneTo(1))) isa Matrix{Int}
    B = similar(Array{Int}, (0:0, 3))
    @test isa(B, OffsetArray{Int, 2})
    @test axes(B) == (0:0, 1:3)

    @test_throws MethodError similar(A, (:,))
    @test_throws MethodError similar(A, (: ,:))
    @test_throws MethodError similar(A, (: ,2))
    @test_throws MethodError similar(A, Float64, (: ,:))
    @test_throws MethodError similar(A, Float64, (: ,2))

    function testsimilar(args...)
        try
           similar(args...)
        catch e
            @test e isa MethodError
            io = IOBuffer()
            showerror(io, e)
            s = split(String(take!(io)),'\n')[1]
            @test occursin(repr(similar), s)
        end
    end

    testsimilar(typeof(A), (:, :))
    testsimilar(typeof(A), (:, 2))
    testsimilar(typeof(A), (:, 1:3))
end

@testset "reshape" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))

    B = reshape(A0, -10:-9, 9:10)
    @test isa(B, OffsetArray{Int,2})
    @test parent(B) === A0
    @test axes(B) == IdentityUnitRange.((-10:-9, 9:10))
    B = reshape(A, -10:-9, 9:10)
    @test isa(B, OffsetArray{Int,2})
    @test pointer(parent(B)) === pointer(A0)
    @test axes(B) == IdentityUnitRange.((-10:-9, 9:10))
    b = reshape(A, -7:-4)
    @test axes(b) == (IdentityUnitRange(-7:-4),)
    @test isa(parent(b), Vector{Int})
    @test pointer(parent(b)) === pointer(parent(A))
    @test parent(b) == A0[:]
    a = OffsetArray(rand(3,3,3), -1:1, 0:2, 3:5)
    # Offset axes are required for reshape(::OffsetArray, ::Val) support
    b = reshape(a, Val(2))
    @test isa(b, OffsetArray{Float64,2})
    @test pointer(parent(b)) === pointer(parent(a))
    @test axes(b) == IdentityUnitRange.((-1:1, 1:9))
    b = reshape(a, Val(4))
    @test isa(b, OffsetArray{Float64,4})
    @test pointer(parent(b)) === pointer(parent(a))
    @test axes(b) == (axes(a)..., IdentityUnitRange(1:1))

    @test reshape(OffsetArray(-1:0, -1:0), :, 1) == reshape(-1:0, 2, 1)
    @test reshape(OffsetArray(-1:2, -1:2), -2:-1, :) == reshape(-1:2, -2:-1, 2)

    @test reshape(OffsetArray(-1:0, -1:0), :) == OffsetArray(-1:0, -1:0)
    @test reshape(A, :) == reshape(A0, :)

    # pop the parent
    B = reshape(A, size(A))
    @test B == A0
    @test parent(B) === A0
    B = reshape(A, (Base.OneTo(2), 2))
    @test B == A0
    @test parent(B) === A0
    B = reshape(A, (2,:))
    @test B == A0
    @test parent(B) === A0

    # julialang/julia #33614
    A = OffsetArray(-1:0, (-2,))
    @test reshape(A, :) === A
    Arsc = reshape(A, :, 1)
    Arss = reshape(A, 2, 1)
    @test Arsc[1,1] == Arss[1,1] == -1
    @test Arsc[2,1] == Arss[2,1] == 0
    @test_throws BoundsError Arsc[0,1]
    @test_throws BoundsError Arss[0,1]
    A = OffsetArray([-1,0], (-2,))
    Arsc = reshape(A, :, 1)
    Arsc[1,1] = 5
    @test first(A) == 5
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
        for s in (IdentityUnitRange(2:3), OffsetArray(2:3, 2:3))
            @test axes(r[s]) == axes(s)
        end
    end
end

@testset "logical indexing" begin
    A0 = [1 3; 2 4]
    A = OffsetArray(A0, (-1,2))

    @test A[A .> 2] == [3,4]
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
end

@testset "reductions" begin
    A = OffsetArray(rand(4,4), (-3,5))
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
        for r in Any[5:100, UnitRange(5.0, 20.0), IdOffsetRange(1:100, 4), IdOffsetRange(4:5), # AbstractUnitRanges
            2:4:14, 1.5:1.0:10.5, # AbstractRanges
            ]

            a = OffsetVector(r, 2);
            @test mapreduce(identity, +, a) == mapreduce(identity, +, r)
            @test mapreduce(x -> x^2, (x,y) -> x, a) == mapreduce(x -> x^2, (x,y) -> x, r)

            b = mapreduce(identity, +, a, dims = 1)
            br = mapreduce(identity, +, r, dims = 1)
            @test no_offset_view(b) == no_offset_view(br)
            @test axes(b, 1) == first(axes(a,1)):first(axes(a,1))

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
                @test axes(b, 1) == first(axes(a,1)):first(axes(a,1))

                b = f(a, dims = 2);
                br = f(r, dims = 2)
                @test no_offset_view(b) == no_offset_view(br)
                @test axes(b, 1) == axes(a,1)
            end

            @test extrema(a) == extrema(r)
        end
    end
end

# v  = OffsetArray([1,1e100,1,-1e100], (-3,))*1000
# v2 = OffsetArray([1,-1e100,1,1e100], (5,))*1000
# @test isa(v, OffsetArray)
# cv  = OffsetArray([1,1e100,1e100,2], (-3,))*1000
# cv2 = OffsetArray([1,-1e100,-1e100,2], (5,))*1000
# @test isequal(cumsum_kbn(v), cv)
# @test isequal(cumsum_kbn(v2), cv2)
# @test isequal(sum_kbn(v), sum_kbn(parent(v)))

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

@testset "fill" begin
    B = fill(5, 1:3, -1:1)
    @test axes(B) == (1:3,-1:1)
    @test all(B.==5)

    B = fill(5, (1:3, -1:1))
    @test axes(B) == (1:3,-1:1)
    @test all(B.==5)

    B = fill(5, 3, -1:1)
    @test axes(B) == (1:3,-1:1)
    @test all(B.==5)
end

@testset "broadcasting" begin
    A = OffsetArray(rand(4,4), (-3,5))

    @test A.+1 == OffsetArray(parent(A).+1, A.offsets)
    @test 2*A == OffsetArray(2*parent(A), A.offsets)
    @test A+A == OffsetArray(parent(A)+parent(A), A.offsets)
    @test A.*A == OffsetArray(parent(A).*parent(A), A.offsets)
end

@testset "@inbounds" begin
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
    OB = OffsetArrays.no_offset_view(B)
    @test axes(OB, 1) == 1:4
    @test collect(OB) == 0:3

    # issue #198
    offax = axes(OffsetVector(1:10, -5), 1)
    noffax = OffsetArrays.no_offset_view(offax)
    @test noffax == -4:5
    @test axes(noffax, 1) == 1:10   # ideally covered by the above, but current it isn't
    @test isa(noffax, AbstractUnitRange)

    r = Base.OneTo(4)
    @test OffsetArrays.no_offset_view(r) isa typeof(r)

    # SubArrays
    A = reshape(1:12, 3, 4)
    V = view(A, OffsetArrays.IdentityUnitRange(2:3), OffsetArrays.IdentityUnitRange(2:3))
    if collect(V) == [5 8; 6 9]   # julia 1.0 has a bug here
        @test OffsetArrays.no_offset_view(V) == [5 8; 6 9]
    end
    V = view(A, OffsetArrays.IdentityUnitRange(2:3), 2)
    @test V != [5;6]
    if collect(V) == [5;6]
        @test OffsetArrays.no_offset_view(V) == [5;6]
    end
    O = OffsetArray(A, -1:1, 0:3)
    V = view(O, 0:1, 1:2)
    @test V == OffsetArrays.no_offset_view(V) == [5 8; 6 9]
    r1, r2 = OffsetArrays.IdOffsetRange(1:3, -2), OffsetArrays.IdentityUnitRange(2:3)
    V = view(O, r1, r2)
    @test V != collect(V)
    @test OffsetArrays.no_offset_view(V) == collect(V)
    V = @view O[:,:]
    @test IndexStyle(A) == IndexStyle(O) == IndexStyle(V) == IndexStyle(OffsetArrays.no_offset_view(V)) == IndexLinear()
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
    @test axes(o, 1) == 0:-1
    @test push!(o, 1) === o
    @test axes(o, 1) == 0:0
    @test o[end] == 1
    @test push!(o, 2, 3) === o
    @test axes(o, 1) == 0:2
    @test o[end-1:end] == [2, 3]
    # pop!
    o = OffsetVector([1, 2, 3], -1)
    @test pop!(o) == 3
    @test axes(o, 1) == 0:1
    # append!
    o = OffsetVector([1, 2, 3], -1)
    append!(o, [4, 5])
    @test axes(o, 1) == 0:4
    # empty!
    o = OffsetVector([1, 2, 3], -1)
    @test empty!(o) === o
    @test axes(o, 1) == 0:-1
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
end

@testset "Adapt" begin
    # We need another storage type, CUDA.jl defines one but we can't use that for CI
    # let's define an appropriate method for SArrays
    Adapt.adapt_storage(::Type{SA}, xs::Array) where SA<:SArray         = convert(SA, xs)   # ambiguity
    Adapt.adapt_storage(::Type{SA}, xs::AbstractArray) where SA<:SArray = convert(SA, xs)
    arr = OffsetArray(rand(3, 3), -1:1, -1:1)
    s_arr = adapt(SMatrix{3,3}, arr)
    @test parent(s_arr) isa SArray
    @test arr == adapt(Array, s_arr)
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
end

include("origin.jl")
