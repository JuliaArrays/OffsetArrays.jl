# For developers

Writing code that supports OffsetArrays is generally fairly straightforward.
The majority of cases can be handled with these tips:

- replace many uses of `size` with `axes`
- replace `1:length(A)` with `eachindex(A)`, or if you need an integer index with `LinearIndices(A)`
- replace explicit allocations like `Array{Int}(undef, size(B))` with `similar(Array{Int}, axes(B))`

More information can be found in [Julia's developer documentation](https://docs.julialang.org/en/v1.0/devdocs/offset-arrays/).
The most subtle issues tend to arise around the axes, and further detail specific to
OffsetArrays.jl follows below.

## Internals

How does OffsetArrays work? The fundamental principle is very simple:
an `OffsetArray` is just a wrapper around a "parent" array, together
with an index offset:

```jldoctest oa; setup=:(using OffsetArrays)
julia> oa = OffsetArray([1 2; 3 4], 0:1, 5:6)
2×2 OffsetArray(::Array{Int64,2}, 0:1, 5:6) with eltype Int64 with indices 0:1×5:6:
 1  2
 3  4

julia> parent(oa)
2×2 Array{Int64,2}:
 1  2
 3  4

julia> oa.offsets
(-1, 4)
```

So `parent(oa)` is the original array we constructed it with, and `oa.offsets` is a tuple,
each entry encoding the index-shift to be applied along the corresponding axis.
When you index `oa[i,j]`, it "translates" the `i,j` indexes back to the parent array's
indexes and then returns the value in the parent.

## The axes of OffsetArrays

```jldoctest oa
julia> axes(oa)
(0:1, 5:6)
```

This looks straightforward, but if you dive deeper you'll notice some complexities:

```jldoctest oa
julia> ax = axes(oa, 2)
5:6

julia> typeof(ax)
OffsetArrays.IdOffsetRange{Int64,Base.OneTo{Int64}}

julia> ax[1]
ERROR: BoundsError: attempt to access 2-element Base.OneTo{Int64} at index [-3]
Stacktrace:
 [1] throw_boundserror(::Base.OneTo{Int64}, ::Int64) at ./abstractarray.jl:538
 [2] getindex at ./range.jl:625 [inlined]
 [3] getindex(::OffsetArrays.IdOffsetRange{Int64,Base.OneTo{Int64}}, ::Int64) at /home/tim/.julia/dev/OffsetArrays/src/axes.jl:139
 [4] top-level scope at none:0

julia> ax[5]
5
```

The axes are of type [`OffsetArrays.IdOffsetRange`](@ref).
`IdOffsetRange`s are essentially OffsetArrays specialized for ranges, with the additional
property that they tend to be their own axes:

```jldoctest oa
julia> ax
5:6

julia> axes(ax)
(5:6,)

julia> axes(ax[ax])
(5:6,)
```

This example of indexing is [idempotent](https://en.wikipedia.org/wiki/Idempotence).
This is a useful characteristic for ensuring the "fundamental axiom" of generalized indexing,
that `a[rng][i] == a[rng[i]]`:

```jldoctest; setup=:(using OffsetArrays)
julia> oa2 = OffsetArray([5, 10, 15, 20], 0:3)
4-element OffsetArray(::Array{Int64,1}, 0:3) with eltype Int64 with indices 0:3:
  5
 10
 15
 20

julia> ax2 = axes(oa2, 1)
0:3

julia> oa2[2]
15

julia> oa2[ax2][2]
15

julia> oa2[ax2[2]]
15
```

`IdOffsetRange`s apply the offset both to the values and the indices of the range, and otherwise preserve the parent range.

!!! warning

    There are circumstances where constructing a specific type of `IdOffsetRange` cannot be supported without changing the axes of the range (see [`OffsetArrays.IdOffsetRange`](@ref).)
    In the future, this package will distinguish between *construction*  and *conversion*:

    - construction (aka, *coercion*) will always succeed, even if it has to change the axes of the result (Examples: `RangeType(rng)`, `typeof(rng1)(rng2)`)
    - conversion will succeed only if it can preserve both the values and the axes (Examples: `convert(RangeType, rng)`, `oftype(rng1, rng2)`)

    While these behave equivalently now (conversion currently performs coercion), developers are encouraged to "future-proof" their code by choosing the behavior appropriate for each usage.
