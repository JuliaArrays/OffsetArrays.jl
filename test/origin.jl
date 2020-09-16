using OffsetArrays: Origin
@testset "Origin" begin
    get_origin(A::AbstractArray) = first.(axes(A))

    @test Origin(0) != Origin((0, ))
    @test Origin(CartesianIndex(1, 2)) === Origin((1, 2)) === Origin(1, 2)

    # 1d
    v = [1, 2]
    @test get_origin(OffsetArray(v, Origin(2))) == (2, )
    ov = OffsetArray(v, -3)
    @test get_origin(OffsetArray(ov, Origin(2))) == (2, )
    @test get_origin(OffsetVector(ov, Origin(2))) == (2, )
    @test get_origin(OffsetArray(ov, Origin((2, )))) == (2, )

    # 2d
    a = [1 2;3 4]
    @test get_origin(OffsetArray(a, Origin(0))) == (0, 0)
    oa = OffsetArray(a, -3, -3)
    @test get_origin(OffsetArray(oa, Origin(0))) == (0, 0)
    @test get_origin(OffsetMatrix(oa, Origin(0))) == (0, 0)
    @test get_origin(OffsetArray(oa, Origin(1, 2))) == (1, 2)

    # 3d
    a = ones(3, 3, 3)
    @test get_origin(OffsetArray(a, Origin(0))) == (0, 0, 0)
    oa = OffsetArray(a, -3, -3, -3)
    @test get_origin(OffsetArray(oa, Origin(0))) == (0, 0, 0)
    @test get_origin(OffsetArray(oa, Origin(1, 2, 3))) == (1, 2, 3)
end
