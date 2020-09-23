"""
    Origin(indices...)
    Origin(origin::Tuple)
    Origin(origin::CartesianIndex)

A helper type to construct OffsetArray with given origin.

The `origin` of an array is defined as the index of its first element, i.e., `first.(axes(A))`.

# Example

```jldoctest; setup=:(using OffsetArrays)
julia> a = [1 2; 3 4];

julia> OffsetArray(a, OffsetArrays.Origin(0, 1))
2×2 OffsetArray(::$(Array{Int,2}), 0:1, 1:2) with eltype $Int with indices 0:1×1:2:
 1  2
 3  4

julia> OffsetArray(a, OffsetArrays.Origin(0)) # short notation for `Origin(0, 0)`
2×2 OffsetArray(::$(Array{Int, 2}), 0:1, 0:1) with eltype $Int with indices 0:1×0:1:
 1  2
 3  4
```
"""
struct Origin{T<:Union{Tuple, Int}}
    index::T
end
Origin(I::NTuple{N, Int}) where N = Origin{typeof(I)}(I)
Origin(I::CartesianIndex) = Origin(I.I)
Origin(I1::Int, In::Int...) = Origin((I1, In...))
# Origin(0) != Origin((0, )) but they work the same with broadcasting
Origin(n::Int) = Origin{Int}(n)

(o::Origin)(A::AbstractArray) = o.index .- first.(axes(A))
