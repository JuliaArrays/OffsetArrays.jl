OffsetArrays.jl
===============

provides the ability to use arbitrary starting indices for arrays in Julia programming language.

```jlcon
julia> using OffsetArrays

julia> y = OffsetArray(Float64, -1:1, -7:7, -128:512, -5:5, -1:1, -3:3, -2:2, -1:1);

julia> y[-1,-7,-128,-5,-1,-3,-2,-1] = 14
14

julia> y[-1,-7,-128,-5,-1,-3,-2,-1] += 5
19.0
```

The main purpose of this package is to simplify translation from Fortran codes intensively using Fortran arrays with negative and zero starting indices,
such as
* the codes accompanying the book _Numerical Solution of Hyperbolic Partial Differential Equations_ by prof. John A. Trangenstein
  Please refer [here](http://www.math.duke.edu/~johnt/) and [here](http://www.cambridge.org/us/academic/subjects/mathematics/differential-and-integral-equations-dynamical-systems-and-co/numerical-solution-hyperbolic-partial-differential-equations)
* Clawpack (stands for *Conservation Laws Package*) by prof. Randall J. LeVeque

  + [classic ClawPack](http://depts.washington.edu/clawpack/)

  + [ClawPack 5.0](http://clawpack.github.io/index.html)

The directory _examples_ contains at the moment a translation of explicit upwind finite difference scheme for scalar law advection equation from
the book _Numerical Solution of Hyperbolic Partial Differential Equations_ by prof. John A. Trangenstein.

+ ./scalar_law/PROGRAM0/main.jl
    The original Fortran arrays  `u(-2:ncells+1), x(0:ncells), flux(0:ncells)` transformed to 1-based Julia arrays by shifting index expressions
```julia
    u    = Array(Float64, ncells+4)
    x    = Array(Float64, ncells+1)
    flux = Array(Float64, ncells+1)
```

+ ./scalar_law/PROGRAM0/main_offset_array.jl
    Exemplary _use case_ of this package
```julia
    u    = OffsetArray(Float64, -2:ncells+1)
    x    = OffsetArray(Float64,  0:ncells)
    flux = OffsetArray(Float64,  0:ncells)
```
