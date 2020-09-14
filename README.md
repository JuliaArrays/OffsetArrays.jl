# OffsetArrays.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaArrays.github.io/OffsetArrays.jl/stable)
[![Build Status](https://travis-ci.org/JuliaArrays/OffsetArrays.jl.svg?branch=master)](https://travis-ci.org/JuliaArrays/OffsetArrays.jl)
[![codecov.io](http://codecov.io/github/JuliaArrays/OffsetArrays.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaArrays/OffsetArrays.jl?branch=master)
[![PkgEval][pkgeval-img]][pkgeval-url]


OffsetArrays provides Julia users with arrays that have arbitrary
indices, similar to those found in some other programming languages
like Fortran.

## Usage

You can construct such arrays as follows:

```julia
OA = OffsetArray(A, axis1, axis2, ...)
```

where you want `OA` to have axes `(axis1, axis2, ...)` and be indexed by values that
fall within these axis ranges. Example:

```julia
using OffsetArrays
julia> A = Float64.(reshape(1:15, 3, 5))
3×5 Matrix{Float64}:
 1.0  4.0  7.0  10.0  13.0
 2.0  5.0  8.0  11.0  14.0
 3.0  6.0  9.0  12.0  15.0

julia> OA = OffsetArray(A, -1:1, 0:4) # OA will have axes (-1:1, 0:4)
3×5 OffsetArray(::Matrix{Float64}, -1:1, 0:4) with eltype Float64 with indices -1:1×0:4:
 1.0  4.0  7.0  10.0  13.0
 2.0  5.0  8.0  11.0  14.0
 3.0  6.0  9.0  12.0  15.0

julia> OA = OffsetArray(A, CartesianIndex(-1, 0):CartesianIndex(1, 4))
3×5 OffsetArray(::Matrix{Float64}, -1:1, 0:4) with eltype Float64 with indices -1:1×0:4:
[...]

julia> OA[-1,0], OA[1,4]
(1.0, 15.0)
```

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/O/OffsetArrays.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html
