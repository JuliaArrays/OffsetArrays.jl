"""
    ro = IdOffsetRange(r::AbstractUnitRange, offset=0)

Construct an "identity offset range". Numerically, `collect(ro) == collect(r) .+ offset`,
with the additional property that `axes(ro, 1) = axes(r, 1) .+ offset`.
When `r` starts at 1, then `ro[i] == i` and even `ro[ro] == ro`,
i.e., it's the "identity," which is the origin of the "Id" in `IdOffsetRange`.

# Examples

The most common case is shifting a range that starts at 1 (either `1:n` or `Base.OneTo(n)`):
```jldoctest; setup=:(import OffsetArrays)
julia> ro = OffsetArrays.IdOffsetRange(1:3, -2)
OffsetArrays.IdOffsetRange(-1:1)

julia> axes(ro, 1)
OffsetArrays.IdOffsetRange(-1:1)

julia> ro[-1]
-1

julia> ro[3]
ERROR: BoundsError: attempt to access 3-element UnitRange{$Int} at index [5]
```

If the range doesn't start at 1, the values may be different from the indices:
```jldoctest; setup=:(import OffsetArrays)
julia> ro = OffsetArrays.IdOffsetRange(11:13, -2)
OffsetArrays.IdOffsetRange(9:11)

julia> axes(ro, 1)     # 11:13 is indexed by 1:3, and the offset is also applied to the axes
OffsetArrays.IdOffsetRange(-1:1)

julia> ro[-1]
9

julia> ro[3]
ERROR: BoundsError: attempt to access 3-element UnitRange{$Int} at index [5]
```

# Extended help

Construction/coercion preserves the (shifted) values of the input range, but may modify
the indexes if required by the specified types. For example,

    r = OffsetArrays.IdOffsetRange{Int,UnitRange{Int}}(3:4)

has `r[1] == 3` and `r[2] == 4`, whereas

    r = OffsetArrays.IdOffsetRange{Int,Base.OneTo{Int}}(3:4)

has `r[3] == 3` and `r[4] == 4`, and `r[1]` would throw a `BoundsError`.
In this latter case, a shift in the axes was needed because `Base.OneTo` ranges
must start with value 1.

!!! warning

    In the future, *conversion* will preserve both the values and
    the indices, throwing an error when this is not achievable. For instance,

        r = convert(OffsetArrays.IdOffsetRange{Int,UnitRange{Int}}, 3:4)

    has `r[1] == 3` and `r[2] == 4` and would satisfy `r == 3:4`, whereas

    ```julia
    julia> convert(OffsetArrays.IdOffsetRange{Int,Base.OneTo{Int}}, 3:4)    # future behavior, not present behavior
    ERROR: ArgumentError: first element must be 1, got 3
    ```

    where the error will arise because the result could not have the same axes as the input.

    An important corollary is that `typeof(r1)(r2)` and `oftype(r1, r2)` will behave differently:
    the first coerces `r2` to be of the type of `r1`, whereas the second converts.
    Developers are urged to future-proof their code by choosing the behavior appropriate for each usage.
"""
struct IdOffsetRange{T<:Integer,I<:AbstractUnitRange{T}} <: AbstractUnitRange{T}
    parent::I
    offset::T

    IdOffsetRange{T,I}(r::I, offset::T) where {T<:Integer,I<:AbstractUnitRange{T}} = new{T,I}(r, offset)
end

# Construction/coercion from arbitrary AbstractUnitRanges
function IdOffsetRange{T,I}(r::AbstractUnitRange, offset::Integer = 0) where {T<:Integer,I<:AbstractUnitRange{T}}
    rc, o = offset_coerce(I, r)
    return IdOffsetRange{T,I}(rc, convert(T, o+offset))
end
function IdOffsetRange{T}(r::AbstractUnitRange, offset::Integer = 0) where T<:Integer
    rc = convert(AbstractUnitRange{T}, r)
    return IdOffsetRange{T,typeof(rc)}(rc, convert(T, offset))
end
IdOffsetRange(r::AbstractUnitRange{T}, offset::Integer = 0) where T<:Integer =
    IdOffsetRange{T,typeof(r)}(r, convert(T, offset))

# Coercion from other IdOffsetRanges
IdOffsetRange{T,I}(r::IdOffsetRange{T,I}) where {T,I} = r
function IdOffsetRange{T,I}(r::IdOffsetRange, offset::Integer = 0) where {T<:Integer,I<:AbstractUnitRange{T}}
    rc, offset_rc = offset_coerce(I, r.parent)
    return IdOffsetRange{T,I}(rc, r.offset + offset + offset_rc)
end
function IdOffsetRange{T}(r::IdOffsetRange, offset::Integer = 0) where T<:Integer
    return IdOffsetRange{T}(r.parent, r.offset + offset)
end
IdOffsetRange(r::IdOffsetRange) = r
IdOffsetRange(r::IdOffsetRange, offset::Integer) = typeof(r)(r.parent, offset + r.offset)

# TODO: uncomment these when Julia is ready
# # Conversion preserves both the values and the indexes, throwing an InexactError if this
# # is not possible.
# Base.convert(::Type{IdOffsetRange{T,I}}, r::IdOffsetRange{T,I}) where {T<:Integer,I<:AbstractUnitRange{T}} = r
# Base.convert(::Type{IdOffsetRange{T,I}}, r::IdOffsetRange) where {T<:Integer,I<:AbstractUnitRange{T}} =
#     IdOffsetRange{T,I}(convert(I, r.parent), r.offset)
# Base.convert(::Type{IdOffsetRange{T,I}}, r::AbstractUnitRange) where {T<:Integer,I<:AbstractUnitRange{T}} =
#     IdOffsetRange{T,I}(convert(I, r), 0)

offset_coerce(::Type{Base.OneTo{T}}, r::Base.OneTo) where T<:Integer = convert(Base.OneTo{T}, r), 0
function offset_coerce(::Type{Base.OneTo{T}}, r::AbstractUnitRange) where T<:Integer
    o = first(r) - 1
    return Base.OneTo{T}(last(r) - o), o
end
# function offset_coerce(::Type{Base.OneTo{T}}, r::IdOffsetRange) where T<:Integer
#     rc, o = offset_coerce(Base.OneTo{T}, r.parent)

# Fallback, specialze this method if `convert(I, r)` doesn't do what you need
offset_coerce(::Type{I}, r::AbstractUnitRange) where I<:AbstractUnitRange =
    convert(I, r), 0

@inline Base.parent(r::IdOffsetRange) = r.parent
@inline Base.axes(r::IdOffsetRange) = (Base.axes1(r),)
@inline Base.axes1(r::IdOffsetRange) = IdOffsetRange(Base.axes1(r.parent), r.offset)
@inline Base.unsafe_indices(r::IdOffsetRange) = (r,)
@inline Base.length(r::IdOffsetRange) = length(r.parent)
Base.reduced_index(i::IdOffsetRange) = typeof(i)(first(i):first(i))
# Workaround for #92 on Julia < 1.4
Base.reduced_index(i::IdentityUnitRange{<:IdOffsetRange}) = typeof(i)(first(i):first(i))
for f in [:firstindex, :lastindex, :first, :last]
    @eval @inline Base.$f(r::IdOffsetRange) = $f(r.parent) + r.offset
end

@inline function Base.iterate(r::IdOffsetRange)
    ret = iterate(r.parent)
    ret === nothing && return nothing
    return (ret[1] + r.offset, ret[2])
end
@inline function Base.iterate(r::IdOffsetRange, i)
    ret = iterate(r.parent, i)
    ret === nothing && return nothing
    return (ret[1] + r.offset, ret[2])
end

@propagate_inbounds Base.getindex(r::IdOffsetRange, i::Integer) = r.parent[i - r.offset] + r.offset
@propagate_inbounds function Base.getindex(r::IdOffsetRange, s::AbstractUnitRange{<:Integer})
    return r.parent[s .- r.offset] .+ r.offset
end
@propagate_inbounds function Base.getindex(r::IdOffsetRange, s::IdentityUnitRange)
    return IdOffsetRange(r.parent[s .- r.offset], r.offset)
end
@propagate_inbounds function Base.getindex(r::IdOffsetRange, s::IdOffsetRange)
    return IdOffsetRange(r.parent[s.parent .+ (s.offset - r.offset)] .+ (r.offset - s.offset), s.offset)
end

# offset-preserve broadcasting
Broadcast.broadcasted(::Base.Broadcast.DefaultArrayStyle{1}, ::typeof(-), r::IdOffsetRange{T}, x::Integer) where T =
    IdOffsetRange{T}(r.parent .- x, r.offset)
Broadcast.broadcasted(::Base.Broadcast.DefaultArrayStyle{1}, ::typeof(+), r::IdOffsetRange{T}, x::Integer) where T =
    IdOffsetRange{T}(r.parent .+ x, r.offset)
Broadcast.broadcasted(::Base.Broadcast.DefaultArrayStyle{1}, ::typeof(+), x::Integer, r::IdOffsetRange{T}) where T =
    IdOffsetRange{T}(x .+ r.parent, r.offset)

Base.show(io::IO, r::IdOffsetRange) = print(io, "OffsetArrays.IdOffsetRange(",first(r), ':', last(r),")")

# Optimizations
@inline Base.checkindex(::Type{Bool}, inds::IdOffsetRange, i::Real) = Base.checkindex(Bool, inds.parent, i - inds.offset)

if VERSION < v"1.5.2"
    # issue 100, 133: IdOffsetRange as another index-preserving case shouldn't comtribute offsets
    # fixed by https://github.com/JuliaLang/julia/pull/37204
    @inline Base.compute_offset1(parent, stride1::Integer, dims::Tuple{Int}, inds::Tuple{IdOffsetRange}, I::Tuple) =
        Base.compute_linindex(parent, I) - stride1*first(inds[1])
end
