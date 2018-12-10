# OffsetArrays.jl


OffsetArrays provides Julia users with arrays that have arbitrary
indices, similar to those found in some other programming languages
like Fortran.

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

## Notes on supporting OffsetArrays

Julia supports generic programming with arrays that doesn't require you to assume that indices start with 1, see the [documentation](http://docs.julialang.org/en/latest/devdocs/offset-arrays/).
