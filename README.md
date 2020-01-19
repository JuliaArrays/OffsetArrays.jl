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
A = reshape(1:15, 3, 5)
println("here is A:")
display(A)
OA = OffsetArray(A, -1:1, 0:4)    # OA will have axes (-1:1, 0:4)
println("here is OA:")
display(OA)
@show OA[-1,0] OA[1,4]
```

which prints out

```
here is A:
3×5 reshape(::UnitRange{Int64}, 3, 5) with eltype Int64:
 1  4  7  10  13
 2  5  8  11  14
 3  6  9  12  15
here is OA:
OffsetArray(reshape(::UnitRange{Int64}, 3, 5), -1:1, 0:4) with eltype Int64 with indices -1:1×0:4:
 1  4  7  10  13
 2  5  8  11  14
 3  6  9  12  15
OA[-1, 0] = 1
OA[1, 4] = 15
```


## Notes on supporting OffsetArrays

There are several "tricks" that make it easier to support arrays with general indexes, see the [documentation](http://docs.julialang.org/en/latest/devdocs/offset-arrays/).


[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/O/OffsetArrays.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html
