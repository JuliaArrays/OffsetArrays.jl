# OffsetArrays.jl

| Documentation (stable) | Documentation (dev) |  Unit Tests | Code coverage | PkgEval |
| :-: | :-: | :-: | :-: | :-: |
| [![][docs-stable-img]][docs-stable-url] |  [![][docs-dev-img]][docs-dev-url] | [![][action-img]][action-url] | [![][codecov-img]][codecov-url] | [![][pkgeval-img]][pkgeval-url] |

## Introduction

OffsetArrays provides Julia users with arrays that have arbitrary
indices, similar to those found in some other programming languages
like Fortran.

An `OffsetArray` is a lightweight wrapper around an `AbstractArray` that shifts its indices.
Generally, indexing into an `OffsetArray` should be as performant as the parent array.

## Usage

There are two ways to construct `OffsetArray`s: by specifying the axes of the array, or
by specifying its origin.

The first way to construct an `OffsetArray` by specifying its axes is:

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

julia> axes(A) # indices of a Matrix start from 1 along each axis
(Base.OneTo(3), Base.OneTo(5))

julia> OA = OffsetArray(A, -1:1, 0:4) # OA will have the axes (-1:1, 0:4)
3×5 OffsetArray(::Matrix{Float64}, -1:1, 0:4) with eltype Float64 with indices -1:1×0:4:
 1.0  4.0  7.0  10.0  13.0
 2.0  5.0  8.0  11.0  14.0
 3.0  6.0  9.0  12.0  15.0

julia> OA[-1,0]
1.0

julia> OA[1,4]
15.0
```

The second way to construct an `OffsetArray` is by specifying the origin, that is, the first index
along each axis. This is particularly useful if one wants, eg., arrays that are 0-indexed as opposed
to 1-indexed.

A convenient way to construct an `OffsetArray` this way is by using `OffsetArrays.Origin`:

```julia
julia> using OffsetArrays: Origin

julia> Origin(0)(A) # indices begin at 0 along all axes
3×5 OffsetArray(::Matrix{Float64}, 0:2, 0:4) with eltype Float64 with indices 0:2×0:4:
 1.0  4.0  7.0  10.0  13.0
 2.0  5.0  8.0  11.0  14.0
 3.0  6.0  9.0  12.0  15.0

julia> Origin(2, 3)(A) # indices begin at 2 along the first axis and 3 along the second
3×5 OffsetArray(::Matrix{Float64}, 2:4, 3:7) with eltype Float64 with indices 2:4×3:7:
 1.0  4.0  7.0  10.0  13.0
 2.0  5.0  8.0  11.0  14.0
 3.0  6.0  9.0  12.0  15.0
```

While the examples here refer to the common case where the parent arrays have indices starting at 1,
this is not necessary. An `OffsetArray` may wrap any array that has integer indices, irrespective of
where the indices begin.

## How to go back to 1-indexed arrays

Certain libraries, such as `LinearAlgebra`, require arrays to be indexed from 1. Passing an `OffsetArray`
with shifted indices would lead to an error here.

```julia
julia> A = Float64.(reshape(1:16, 4, 4));

julia> AO = Origin(0)(A);

julia> using LinearAlgebra

julia> Diagonal(AO)
ERROR: ArgumentError: offset arrays are not supported but got an array with index other than 1
```

The way to obtain a `1`-indexed array from an `OffsetArray` is by using `OffsetArrays.no_offset_view`.

An example of this is:

```julia
julia> OffsetArrays.no_offset_view(AO)
4×4 Matrix{Float64}:
 1.0  5.0   9.0  13.0
 2.0  6.0  10.0  14.0
 3.0  7.0  11.0  15.0
 4.0  8.0  12.0  16.0
```

This may now be passed to `LinearAlgebra`:

```julia
julia> D = Diagonal(OffsetArrays.no_offset_view(AO))
4×4 Diagonal{Float64, Vector{Float64}}:
 1.0   ⋅     ⋅     ⋅
  ⋅   6.0    ⋅     ⋅
  ⋅    ⋅   11.0    ⋅
  ⋅    ⋅     ⋅   16.0
```

If we want to restore the original indices of `AO`, we may wrap an `OffsetArray` around the `Diagonal` as:

```julia
julia> Origin(AO)(D)
4×4 OffsetArray(::Diagonal{Float64, Vector{Float64}}, 0:3, 0:3) with eltype Float64 with indices 0:3×0:3:
 1.0   ⋅     ⋅     ⋅
  ⋅   6.0    ⋅     ⋅
  ⋅    ⋅   11.0    ⋅
  ⋅    ⋅     ⋅   16.0
```

Here, `Origin(AO)` is able to automatically infer and use the indices of `AO`.

<!-- badges -->

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/O/OffsetArrays.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html

[action-img]: https://github.com/JuliaArrays/OffsetArrays.jl/workflows/Unit%20test/badge.svg
[action-url]: https://github.com/JuliaArrays/OffsetArrays.jl/actions

[codecov-img]: https://codecov.io/github/JuliaArrays/OffsetArrays.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/gh/JuliaArrays/OffsetArrays.jl

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://juliaarrays.github.io/OffsetArrays.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://juliaarrays.github.io/OffsetArrays.jl/dev/
