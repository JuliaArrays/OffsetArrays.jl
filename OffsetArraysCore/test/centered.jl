module CenteredTests

using OffsetArraysCore
using OffsetArraysCore: center, centered
using Test

@testset "center/centered" begin
    @testset "center" begin
        A = reshape(collect(1:9), 3, 3)
        c = center(A)
        @test c == (2, 2)
        @test A[c...] == 5
        @test center(A, RoundDown) == center(A, RoundUp)

        A = reshape(collect(1:6), 2, 3)
        c = center(A)
        @test center(A, RoundDown) == c
        @test c == (1, 2)
        @test A[c...] == 3
        c = center(A, RoundUp)
        @test c == (2, 2)
        @test A[c...] == 4
    end

    @testset "centered" begin
        A = reshape(collect(1:9), 3, 3)
        Ao = centered(A)
        @test centered(Ao) === Ao
        @test centered(Ao, center(Ao)) === Ao
        @test typeof(Ao) <: OffsetArray
        @test parent(Ao) === A
        @test Ao.offsets == (-2, -2)
        @test Ao[0, 0] == 5

        A = reshape(collect(1:6), 2, 3)
        Ao = centered(A)
        @test centered(A, center(A, RoundDown)) == Ao
        @test typeof(Ao) <: OffsetArray
        @test parent(Ao) === A
        @test Ao.offsets == (-1, -2)
        @test Ao[0, 0] == 3
        Ao = centered(A, center(A, RoundUp))
        @test typeof(Ao) <: OffsetArray
        @test parent(Ao) === A
        @test Ao.offsets == (-2, -2)
        @test Ao[0, 0] == 4

        A = reshape(collect(1:9), 3, 3)
        Ao = OffsetArray(A, -1, -1)
        Aoo = centered(Ao)
        @test parent(Aoo) === A # there will be only one OffsetArray wrapper
        @test Aoo.offsets == (-2, -2)
        @test Aoo[0, 0] == 5

        A = reshape(collect(1:9), 3, 3)
        Aoo = centered(A, CartesianIndex(2,2))
        c = (0,0)
        i = CartesianIndex(c...)
        @test Aoo[i] == Aoo[c...]

    end
end

@info "Following deprecations are expected"
@testset "deprecations" begin
    A = reshape(collect(1:9), 3, 3)
    @test centered(A, RoundDown) == centered(A, RoundUp)
end

end