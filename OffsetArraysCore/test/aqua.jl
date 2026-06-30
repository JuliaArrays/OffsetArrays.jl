module AquaTests

import Aqua
import OffsetArraysCore
using Test

@testset "Project meta quality checks" begin
    Aqua.test_all(OffsetArraysCore, piracies=false)
end

end