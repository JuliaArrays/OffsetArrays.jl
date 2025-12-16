module DocStringTests

using OffsetArraysCore
using Test
using Documenter

DocMeta.setdocmeta!(OffsetArraysCore, :DocTestSetup, :(using OffsetArraysCore); recursive=true)

@testset "Project meta quality checks" begin
    if VERSION >= v"1.2"
        doctest(OffsetArraysCore, manual = false)
    end
end

end