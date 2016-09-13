# OffsetArrays.jl


OffsetArrays provides Julia users with arrays that have arbitrary
indices, similar to those found in some other programming languages
like Fortran.

```julia
julia> using OffsetArrays

julia> y = OffsetArray(Float64, -1:1, -7:7, -128:512, -5:5, -1:1, -3:3, -2:2, -1:1);

julia> summary(y)
"OffsetArrays.OffsetArray{Float64,8,Array{Float64,8}} with indices -1:1×-7:7×-128:512×-5:5×-1:1×-3:3×-2:2×-1:1"

julia> y[-1,-7,-128,-5,-1,-3,-2,-1] = 14
14

julia> y[-1,-7,-128,-5,-1,-3,-2,-1] += 5
19.0
```

Support for such arrays is based on new functionality in Julia v0.5,
and the modern package is not usable on earlier Julia versions.

## Special note for Julia 0.5

During the transition during Julia 0.5 to allowing arrays with
arbitrary indices, certain operations (like `size` and `length`) are
deliberately unsupported; see the
[documentation](http://docs.julialang.org/en/latest/devdocs/offset-arrays/)
describing the reasoning behind this decision. The general
recommendation is to use `indices` and `linearindices` for most
operations where `size` and `length` would have formerly been used.

If your package makes use of OffsetArrays, you can also add the
following internal convenience definitions:

```jl
_size(A::AbstractArray) = map(length, indices(A))
_size(A) = size(A)

_length(A::AbstractArray) = length(linearindices(A))
_length(A) = length(A)
```

These versions should work for all types.

## Performance optimization with @unsafe

Also during Julia 0.5, `@inbounds` will not remove the internal
bounds-checking that happens when using an OffsetArray. Until this
changes, you can often accomplish the same task with `@unsafe`:

```jl
v = OffsetArray(rand(1000), 0:999)
@unsafe for i in indices(v, 1)
    v[i] += 1
end
```

With such annotation, `OffsetArray`s are as performant as `Array`s.

`@unsafe` is not as powerful as `@inbounds`, and it is possible to set
up situations where it fails to remove bounds checks.  However, it
suffices for many uses.

# Comparison with Fortran

The original purpose of the package was to simplify translation of Fortran codes, say
* the codes accompanying the book _Numerical Solution of Hyperbolic Partial Differential Equations_ by prof. John A. Trangenstein
  Please refer [here](http://www.math.duke.edu/~johnt/) and [here](http://www.cambridge.org/us/academic/subjects/mathematics/differential-and-integral-equations-dynamical-systems-and-co/numerical-solution-hyperbolic-partial-differential-equations)
* Clawpack (stands for *Conservation Laws Package*) by prof. Randall J. LeVeque

  + [classic ClawPack](http://depts.washington.edu/clawpack/)

  + [ClawPack 5.0](http://clawpack.github.io/index.html)

The directory _examples_ contains at the moment a translation of explicit upwind finite difference scheme for scalar law advection equation from
the book _Numerical Solution of Hyperbolic Partial Differential Equations_ by prof. John A. Trangenstein.

+ examples/scalar_law/PROGRAM0/main.jl
    The original Fortran arrays  `u(-2:ncells+1), x(0:ncells), flux(0:ncells)` transformed to 1-based Julia arrays by shifting index expressions
```julia
    u    = Array(Float64, ncells+4)
    x    = Array(Float64, ncells+1)
    flux = Array(Float64, ncells+1)
```

+ examples/scalar_law/PROGRAM0/main_offset_array.jl
    Exemplary _use case_ of this package
```julia
    u    = OffsetArray(Float64, -2:ncells+1)
    x    = OffsetArray(Float64,  0:ncells)
    flux = OffsetArray(Float64,  0:ncells)
```

+ UPDATE:
    Added
    + examples/scalar_law/PROGRAM0/main_sub.jl

    see more details here
    + [room for performance improvement for SubArray #5117](https://github.com/JuliaLang/julia/issues/5117)

+ Timing for baseline with `@unsafe` annotation:
```sh
  0.103402 seconds (21 allocations: 235.531 KB)
```

+ Timing for `OffsetArray` version with `@unsafe` annotation:
```sh
  0.103987 seconds (16 allocations: 235.094 KB)
```

+ Timing for `sub` macro version with `@unsafe` annotation (doesn't work without `@unsafe`):
```sh
  0.105967 seconds (217 allocations: 268.094 KB)
```
Only the 2nd timing after warming up is given.
