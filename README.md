# OffsetArrays.jl

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

OffsetArrays works for arbitrary dimensionality:

```julia
julia> using OffsetArrays

julia> y = OffsetArray{Float64}(undef, -1:1, -7:7, -128:512, -5:5, -1:1, -3:3, -2:2, -1:1);

julia> summary(y)
"OffsetArrays.OffsetArray{Float64,8,Array{Float64,8}} with indices -1:1×-7:7×-128:512×-5:5×-1:1×-3:3×-2:2×-1:1"

julia> y[-1,-7,-128,-5,-1,-3,-2,-1] = 14
14

julia> y[-1,-7,-128,-5,-1,-3,-2,-1] += 5
19.0
```

## Example: Relativistic Notation
Suppose we have a position vector `r = [:x, :y, :z]` which is naturally one-based, ie. `r[1] == :x`, `r[2] == :y`,  `r[3] == :z` and we also want to construct a relativistic position vector which includes time as the 0th component. This can be done with OffsetArrays like

```julia
julia> using OffsetArrays

julia> r = [:x, :y, :z];

julia> x = OffsetVector([:t, r...], 0:3)
OffsetArray(::Array{Symbol,1}, 0:3) with eltype Symbol with indices 0:3:
 :t
 :x
 :y
 :z

julia> x[0]
:t

julia> x[1:3]
3-element Array{Symbol,1}:
 :x
 :y
 :z
```

## Example: Polynomials
Suppose one wants to represent the Laurent polynomial
```
6/x + 5 - 2*x + 3*x^2 + x^3
```
in julia. The coefficients of this polynomial are a naturally `-1` based list, since the `n`th element of the list
(counting from `-1`) `6, 5, -2, 3, 1` is the coefficient corresponding to the `n`th power of `x`. This Laurent polynomial can be evaluated at say `x = 2` as follows.
```julia
julia> using OffsetArrays

julia> coeffs = OffsetVector([6, 5, -2, 3, 1], -1:3)
OffsetArray(::Array{Int64,1}, -1:3) with eltype Int64 with indices -1:3:
  6
  5
 -2
  3
  1

julia> polynomial(x, coeffs) = sum(coeffs[n]*x^n for n in eachindex(coeffs))
polynomial (generic function with 1 method)

julia> polynomial(2.0, coeffs)
24.0
```
Notice our use of the `eachindex` function which does not assume that the given array starts at `1`.

## Notes on supporting OffsetArrays

There are several "tricks" that make it easier to support arrays with general indexes, see the [documentation](http://docs.julialang.org/en/latest/devdocs/offset-arrays/).


[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/O/OffsetArrays.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html
