# Comparison with Fortran

The original purpose of the package was to simplify translation of Fortran codes, say
* the codes accompanying the book _Numerical Solution of Hyperbolic Partial Differential Equations_ by prof. John A. Trangenstein
  Please refer [here](http://www.math.duke.edu/~johnt/) and [here](http://www.cambridge.org/us/academic/subjects/mathematics/differential-and-integral-equations-dynamical-systems-and-co/numerical-solution-hyperbolic-partial-differential-equations)
* Clawpack (stands for *Conservation Laws Package*) by prof. Randall J. LeVeque

  + [classic ClawPack](http://depts.washington.edu/clawpack/)

  + [ClawPack 5.0](http://clawpack.github.io/index.html)

see also [translation to Julia of Fortran codes from ClawPack](https://github.com/alsam/Claw.jl)

The directory _examples_ of [hyperbolic_PDE](https://github.com/alsam/hyperbolic_PDE.jl) contains at the moment a translation of explicit upwind finite difference scheme for scalar law advection equation from
the book _Numerical Solution of Hyperbolic Partial Differential Equations_ by prof. John A. Trangenstein.

+ examples/scalar_law/PROGRAM0/main.jl
    The original Fortran arrays  `u(-2:ncells+1), x(0:ncells), flux(0:ncells)` transformed to 1-based Julia arrays by shifting index expressions
```julia
    u    = Array{Float64}(undef, ncells+4)
    x    = Array{Float64}(undef, ncells+1)
    flux = Array{Float64}(undef, ncells+1)
```{+ examp}(undef,es/scalar_law/PROGRAM0/main_offset_array.jl
    Exemplary _use case_ of this package
```julia
    u    = OffsetArray{Float64}(undef, -2:ncells+1)
    x    = OffsetArray{Float64}(undef,  0:ncells)
    flux = OffsetArray{Float64}(undef,  0:ncells)
```

+ Timings for baseline with `@inbounds` annotation:
```sh
  0.103402 seconds (21 allocations: 235.531 KB)
```

+ Timing for `OffsetArray` version with `@inbounds` annotation:
```sh
  0.103987 seconds (16 allocations: 235.094 KB)
```

+ UPDATE:
    Added
    + examples/scalar_law/PROGRAM1/...
```sh
[~/w/m/O/e/s/PROGRAM1] $ julia linaddmain.jl --cells=10000 --runs=3                                                                           ms  master|✚ 1…
  0.672295 seconds (42.90 k allocations: 1.990 MB)
  0.509693 seconds (18 allocations: 313.281 KB)
  0.512243 seconds (18 allocations: 313.281 KB)
[~/w/m/O/e/s/PROGRAM1] $ julia linaddmain.jl --cells=100000 --runs=3                                                                      6134ms  master|✚ 1…
  7.270463 seconds (42.90 k allocations: 4.736 MB)
  7.177485 seconds (18 allocations: 3.053 MB)
  7.248687 seconds (18 allocations: 3.053 MB)
```

    Fortran timing for `100000` number of cells
```sh
real    0m4.781s
user    0m4.760s
sys     0m0.000s

```

    That is 7.2s for Julia script vs. 4.8s for Fortran.
