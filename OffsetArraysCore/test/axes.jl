module AxesTests

using OffsetArraysCore
using OffsetArraysCore: IdOffsetRange, IdentityUnitRange, no_offset_view
using Test

function same_value(r1, r2)
    length(r1) == length(r2) || return false
    for (v1, v2) in zip(r1, r2)
        v1 == v2 || return false
    end
    return true
end

no_offset_axes(x, d) = no_offset_view(axes(x, d))
no_offset_axes(x) = map(no_offset_view, axes(x))

@testset "IdOffsetRange" begin

    function check_indexed_by(r, rindx)
        for i in rindx
            r[i]
        end
        @test_throws BoundsError r[minimum(rindx)-1]
        @test_throws BoundsError r[maximum(rindx)+1]
        return nothing
    end

    ro = IdOffsetRange(Base.OneTo(3))
    rs = IdOffsetRange(3:5, -2)
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
    @test @inferred(IdOffsetRange{Int}(ro))   === ro
    @test @inferred(IdOffsetRange{Int16}(ro)) === IdOffsetRange(Base.OneTo(Int16(3)))
    @test @inferred(IdOffsetRange(ro))        === ro
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
    r = IdOffsetRange{Int,UnitRange{Int}}(3:4)
    @test same_value(r, 3:4)
    check_indexed_by(r, 1:2)
    r = IdOffsetRange{Int,Base.OneTo{Int}}(3:4)
    @test same_value(r, 3:4)
    check_indexed_by(r, 3:4)
    r = IdOffsetRange{Int,Base.OneTo{Int}}(3:4, -2)
    @test same_value(r, 1:2)
    check_indexed_by(r, 1:2)

    r = IdOffsetRange{Int32, Base.OneTo{Int32}}(Base.OneTo(Int64(2)), 3)
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

    # eltype coercion through the AbstractUnitRange constructor
    ro = IdOffsetRange(Base.OneTo(3))
    @test @inferred(AbstractUnitRange{Int}(ro)) === ro
    rb = IdOffsetRange(Base.OneTo(big(3)))
    @test @inferred(AbstractUnitRange{Int}(rb)) === IdOffsetRange(Base.OneTo(3))

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

    p = IdOffsetRange(1:3, 2)
    q = IdOffsetRange(values = p .- 2, indices = p)
    @test same_value(q, 1:3)
    check_indexed_by(q, p)

    @testset for indices in Any[Base.OneTo(3), IdentityUnitRange(Base.OneTo(3))]
        p = IdOffsetRange(values = IdOffsetRange(1:3, 2), indices = indices)
        @test same_value(p, 3:5)
        check_indexed_by(p, 1:3)
        q = IdOffsetRange(values = Base.OneTo(3), indices = indices)
        @test same_value(q, 1:3)
        @test q isa IdOffsetRange{Int, Base.OneTo{Int}}
    end

    # conversion preserves both the values and the axes, throwing an error if this is not possible
    @test @inferred(oftype(ro, ro)) === ro
    @test @inferred(convert(IdOffsetRange{Int}, ro)) === ro
    @test @inferred(convert(IdOffsetRange{Int}, rs)) === rs
    @test @inferred(convert(IdOffsetRange{Int16}, ro)) === IdOffsetRange(Base.OneTo(Int16(3)))
    r2 = @inferred(oftype(rs, ro))
    @test typeof(r2) === typeof(rs)
    @test same_value(r2, 1:3)
    check_indexed_by(r2, 1:3)
    # These two broken tests can be fixed by uncommenting the `convert` definitions
    # in axes.jl, but unfortunately Julia may not quite be ready for this. (E.g. `reinterpretarray.jl`)
    @test_broken try oftype(ro, rs); false catch err true end  # replace with line below
    # @test_throws ArgumentError oftype(ro, rs)
    @test @inferred(oftype(ro, Base.OneTo(2))) === IdOffsetRange(Base.OneTo(2))
    @test @inferred(oftype(ro, 1:2)) === IdOffsetRange(Base.OneTo(2))
    @test_broken try oftype(ro, 3:4); false catch err true end
    # @test_throws ArgumentError oftype(ro, 3:4)

    # broadcasting behavior with scalars (issue #104)
    r3 = (1 .+ IdOffsetRange(3:5, -1) .+ 1) .- 1
    @test r3 isa IdOffsetRange
    @test same_value(r3, 3:5)
    check_indexed_by(r3, axes(r3,1))

    r = IdOffsetRange(3:5, -1)
    rc = copyto!(similar(r), r)
    n = big(typemax(Int))
    @test @inferred(broadcast(+, r, n)) == @inferred(broadcast(+, n, r)) == rc .+ n
    @test @inferred(broadcast(-, r, n)) == rc .- n
    @test @inferred(broadcast(big, r)) == big.(rc)
    for n in Any[2, big(typemax(Int))]
        @test @inferred(broadcast(+, r, n)) == @inferred(broadcast(+, n, r)) == rc .+ n
    end

    @testset "Idempotent indexing" begin
        @testset "Indexing into an IdOffsetRange" begin
            r = IdOffsetRange(3:5, -1)
            # Indexing with IdentityUnitRange
            s = IdentityUnitRange(0:2)
            @test axes(r[s]) == axes(s)
            for i in eachindex(s)
                @test r[s[i]] == r[s][i]
            end

            # Indexing with IdOffsetRange
            s = IdOffsetRange(-4:-2, 4)
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
            r = IdOffsetRange(3:5, -1)
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

    @testset "reduced_indices" begin
        a = reshape(1:24, 2, 3, 4)
        sa = OffsetArray(a, (2, 3, 4));
        @testset for dim in 1:ndims(sa)
            sasum = sum(sa, dims = dim)
            @test parent(sasum) == sum(a, dims = dim)
            find = firstindex(sa, dim)
            @test no_offset_axes(sasum, dim) == find:find
        end
    end

    @testset "conversion to AbstractUnitRange" begin
        r = IdOffsetRange(1:2)
        @test AbstractUnitRange{Int}(r) === r
        r2 = IdOffsetRange(big(1):big(2))
        @test AbstractUnitRange{Int}(r2) === r
        @test AbstractUnitRange{BigInt}(r2) === r2

        if v"1.5" < VERSION
            @test OrdinalRange{Int,Int}(r2) === r
            @test OrdinalRange{BigInt,BigInt}(r2) === r2
        end
    end

    @testset "Bool IdOffsetRange (issue #223)" begin
        for b1 in [false, true], b2 in [false, true]
            r = IdOffsetRange(b1:b2)
            @test first(r) === b1
            @test last(r) === b2
        end
        @test_throws ArgumentError IdOffsetRange(true:true, true)
        @test_throws ArgumentError IdOffsetRange{Bool,UnitRange{Bool}}(true:true, true)
        @test_throws ArgumentError IdOffsetRange{Bool,IdOffsetRange{Bool,UnitRange{Bool}}}(IdOffsetRange(true:true), true)
    end

    @testset "Logical indexing" begin
        @testset "indexing with a single bool" begin
            r = IdOffsetRange(1:2)
            @test_throws ArgumentError r[true]
            @test_throws ArgumentError r[false]
        end
        @testset "indexing with a Bool UnitRange" begin
            r = IdOffsetRange(1:0)

            @test r[true:false] == 1:0
            @test r[true:false] == collect(r)[true:false]
            @test_throws BoundsError r[true:true]
            @test_throws BoundsError r[false:false]
            @test_throws BoundsError r[false:true]

            r = IdOffsetRange(1:1)

            @test r[true:true] == 1:1
            @test r[true:true] == collect(r)[true:true]

            @test r[false:false] == 1:0
            @test r[false:false] == collect(r)[false:false]

            @test_throws BoundsError r[true:false]
            @test_throws BoundsError r[false:true]

            r = IdOffsetRange(1:2)

            @test r[false:true] == 2:2
            @test r[false:true] == collect(r)[false:true]

            @test_throws BoundsError r[true:true]
            @test_throws BoundsError r[true:false]
            @test_throws BoundsError r[false:false]
        end
        @testset "indexing with a Bool IdOffsetRange" begin
            # bounds-checking requires the axes of the indices to match that of the array
            function testlogicalindexing(r, r2)
                r3 = r[r2];
                @test no_offset_view(r3) == collect(r)[collect(r2)]
            end

            r = IdOffsetRange(10:9)
            r2 = IdOffsetRange(true:false)
            testlogicalindexing(r, r2)

            r = IdOffsetRange(10:10)
            r2 = IdOffsetRange(false:false)
            testlogicalindexing(r, r2)
            r2 = IdOffsetRange(true:true)
            testlogicalindexing(r, r2)

            r = IdOffsetRange(10:10, 1)
            r2 = IdOffsetRange(false:false, 1) # effectively true:true with indices 2:2
            testlogicalindexing(r, r2)

            r = IdOffsetRange(10:11)
            r2 = IdOffsetRange(false:true)
            testlogicalindexing(r, r2)
        end
        @testset "indexing with a Bool StepRange" begin
            r = IdOffsetRange(1:0)

            @test r[true:true:false] == 1:1:0
            @test_throws BoundsError r[true:true:true]
            @test_throws BoundsError r[false:true:false]
            @test_throws BoundsError r[false:true:true]

            r = IdOffsetRange(1:1)

            @test r[true:true:true] == 1:1:1
            @test r[true:true:true] == collect(r)[true:true:true]
            @test axes(r[true:true:true], 1) == 1:1

            @test r[false:true:false] == 1:1:0
            @test r[false:true:false] == collect(r)[false:true:false]

            # StepRange{Bool,Int}
            s = StepRange(true, 1, true)
            @test r[s] == 1:1:1
            @test r[s] == collect(r)[s]

            s = StepRange(true, 2, true)
            @test r[s] == 1:1:1
            @test r[s] == collect(r)[s]

            s = StepRange(false, 1, false)
            @test r[s] == 1:1:0
            @test r[s] == collect(r)[s]

            s = StepRange(false, 2, false)
            @test r[s] == 1:1:0
            @test r[s] == collect(r)[s]

            @test_throws BoundsError r[true:true:false]
            @test_throws BoundsError r[false:true:true]

            r = IdOffsetRange(1:2)

            @test r[false:true:true] == 2:1:2
            @test r[false:true:true] == collect(r)[false:true:true]

            # StepRange{Bool,Int}
            s = StepRange(false, 1, true)
            @test r[s] == 2:1:2
            @test r[s] == collect(r)[s]

            @test_throws BoundsError r[true:true:true]
            @test_throws BoundsError r[true:true:false]
            @test_throws BoundsError r[false:true:false]
        end
    end

    @testset "iteration" begin
        # parent has Base.OneTo axes
        A = OffsetArray(ones(7), 4:10)
        ax = axes(A, 1)
        ind, st = iterate(ax)
        @test A[ind] == A[4]
        ind, st = iterate(ax, st)
        @test A[ind] == A[5]

        # parent doesn't have Base.OneTo axes
        B = @view A[:]
        C = OffsetArray(B, 0)
        ax = axes(C, 1)
        ind, st = iterate(ax)
        @test C[ind] == C[4]
        ind, st = iterate(ax, st)
        @test C[ind] == C[5]
    end
end

end