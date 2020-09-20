@testset "IdOffsetRange" begin
    function same_value(r1, r2)
        length(r1) == length(r2) || return false
        for (v1, v2) in zip(r1, r2)
            v1 == v2 || return false
        end
        return true
    end
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
end
