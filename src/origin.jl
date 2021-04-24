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
struct Origin{T<:Union{Tuple{Vararg{Int}}, Int}}
    index::T
end
Origin(I::Tuple{Vararg{Int}}) = Origin{typeof(I)}(I)
Origin(I::Tuple{Vararg{Number}}) = Origin(map(Int, I))
Origin(I::CartesianIndex) = Origin(Tuple(I))
Origin(I::Number...) = Origin(I)
# Origin(0) != Origin((0, )) but they work the same with broadcasting
Origin(n::Number) = Origin{Int}(Int(n))

(o::Origin)(A::AbstractArray) = o.index .- first.(axes(A))

Base.Broadcast.broadcastable(o::Origin) = Ref(o)
