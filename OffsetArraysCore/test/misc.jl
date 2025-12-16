module MiscTests

using OffsetArraysCore
using Test

@testset "misc" begin
    @test OffsetArraysCore._subtractoffset(Base.OneTo(2), 1) isa AbstractUnitRange{Int}
    @test OffsetArraysCore._subtractoffset(Base.OneTo(2), 1) == 0:1
    @test OffsetArraysCore._subtractoffset(3:2:9, 1) isa AbstractRange{Int}
    @test OffsetArraysCore._subtractoffset(3:2:9, 1) == 2:2:8

    @test OffsetArraysCore._addoffset(Base.OneTo(2), 1) isa AbstractUnitRange{Int}
    @test OffsetArraysCore._addoffset(Base.OneTo(2), 1) == 2:3
    @test OffsetArraysCore._addoffset(3:2:9, 1) isa AbstractRange{Int}
    @test OffsetArraysCore._addoffset(3:2:9, 1) == 4:2:10
end

end