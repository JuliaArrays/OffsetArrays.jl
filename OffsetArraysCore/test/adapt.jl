module AdaptTests

using OffsetArraysCore
using Adapt
using StaticArrays
using Test

@testset "Adapt" begin
    # We need another storage type, CUDA.jl defines one but we can't use that for CI
    # let's define an appropriate method for SArrays
    Adapt.adapt_storage(::Type{SA}, xs::Array) where SA<:SArray         = convert(SA, xs)   # ambiguity
    Adapt.adapt_storage(::Type{SA}, xs::AbstractArray) where SA<:SArray = convert(SA, xs)
    arr = OffsetArray(rand(3, 3), -1:1, -1:1)
    s_arr = adapt(SMatrix{3,3}, arr)
    @test parent(s_arr) isa SArray
    @test arr == adapt(Array, s_arr)

    arr2 = OffsetArray(view(rand(5, 5), 2:4, 2:4), -1:1, -1:1)

    if isdefined(Adapt, :parent_type)
        @test Adapt.parent_type(typeof(arr)) == typeof(arr.parent)
        @test Adapt.unwrap_type(typeof(arr)) == typeof(arr.parent)
        @test Adapt.unwrap_type(typeof(arr2)) == typeof(arr.parent)
    end
end

end