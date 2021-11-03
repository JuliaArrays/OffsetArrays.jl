module OffsetArrays

using Base: tail, @propagate_inbounds
@static if !isdefined(Base, :IdentityUnitRange)
    const IdentityUnitRange = Base.Slice
else
    using Base: IdentityUnitRange
end

export OffsetArray, OffsetMatrix, OffsetVector

const IIUR = IdentityUnitRange{<:AbstractUnitRange{<:Integer}}

include("axes.jl")
include("utils.jl")
include("origin.jl")

# Technically we know the length of CartesianIndices but we need to convert it first, so here we
# don't put it in OffsetAxisKnownLength.
const OffsetAxisKnownLength = Union{Integer, AbstractUnitRange}
const OffsetAxis = Union{OffsetAxisKnownLength, Colon}
const ArrayInitializer = Union{UndefInitializer, Missing, Nothing}

## OffsetArray
"""
    OffsetArray(A, indices...)

Return an `AbstractArray` that shares element type and size with the first argument but
uses the supplied `indices` to infer its axes. If all the indices are `AbstractUnitRange`s then
these are directly used as the axis span along each dimension. Refer to the examples below for other
permissible types.

Alternatively it's possible to specify the coordinates of one corner of the array
and have the axes be computed automatically from the size of `A`.
This constructor makes it convenient to shift to
an arbitrary starting index along each axis, for example to a zero-based indexing scheme followed by
arrays in languages such as C and Python.
See [`Origin`](@ref) and the examples below for this usage.

# Example: offsets

There are two types of `indices`: integers and ranges-like types.

Integers are recognized as offsets, where `0` means no offsets are applied:

```jldoctest; setup=:(using OffsetArrays)
julia> A = OffsetArray(reshape(1:6, 2, 3), -1, -2)
2×3 OffsetArray(reshape(::UnitRange{$Int}, 2, 3), 0:1, -1:1) with eltype $Int with indices 0:1×-1:1:
 1  3  5
 2  4  6

julia> A[0, 1]
5
```

Examples of range-like types are: `UnitRange` (e.g, `-1:2`), `CartesianIndices`,
and `Colon()` (or concisely `:`). A `UnitRange` specifies the axis span along one particular dimension,
`CartesianIndices` specify the axis spans along multiple dimensions, and a `Colon` is a placeholder
that specifies that the `OffsetArray` shares its axis with its parent along that dimension.

```jldoctest; setup=:(using OffsetArrays)
julia> OffsetArray(reshape(1:6, 2, 3), 0:1, -1:1)
2×3 OffsetArray(reshape(::UnitRange{$Int}, 2, 3), 0:1, -1:1) with eltype $Int with indices 0:1×-1:1:
 1  3  5
 2  4  6

julia> OffsetArray(reshape(1:6, 2, 3), :, -1:1) # : as a placeholder to indicate that no offset is to be applied to the first dimension
2×3 OffsetArray(reshape(::UnitRange{$Int}, 2, 3), 1:2, -1:1) with eltype $Int with indices 1:2×-1:1:
 1  3  5
 2  4  6
```

Use `CartesianIndices` to specify the coordinates of two diagonally opposite corners:

```jldoctest; setup=:(using OffsetArrays)
julia> OffsetArray(reshape(1:6, 2, 3), CartesianIndex(0, -1):CartesianIndex(1, 1))
2×3 OffsetArray(reshape(::UnitRange{$Int}, 2, 3), 0:1, -1:1) with eltype $Int with indices 0:1×-1:1:
 1  3  5
 2  4  6
```

Integers and range-like types may not be combined in the same call:

```julia
julia> OffsetArray(reshape(1:6, 2, 3), 0, -1:1)
ERROR: [...]
```

# Example: origin

[`OffsetArrays.Origin`](@ref) may be used to specify the origin of the OffsetArray. The term origin here
refers to the corner with the lowest values of coordinates, such as the left edge for an `AbstractVector`,
the bottom left corner for an `AbstractMatrix` and so on. The coordinates of the origin sets the starting
index of the array along each dimension.

```jldoctest; setup=:(using OffsetArrays)
julia> a = [1 2; 3 4];

julia> OffsetArray(a, OffsetArrays.Origin(0, 1))
2×2 OffsetArray(::$(Array{Int,2}), 0:1, 1:2) with eltype $Int with indices 0:1×1:2:
 1  2
 3  4

julia> OffsetArray(a, OffsetArrays.Origin(0)) # set the origin to zero along each dimension
2×2 OffsetArray(::$(Array{Int, 2}), 0:1, 0:1) with eltype $Int with indices 0:1×0:1:
 1  2
 3  4
```


"""
struct OffsetArray{T,N,AA<:AbstractArray{T,N}} <: AbstractArray{T,N}
    parent::AA
    offsets::NTuple{N,Int}
    @inline function OffsetArray{T, N, AA}(parent::AA, offsets::NTuple{N, Int}; checkoverflow = true) where {T, N, AA<:AbstractArray{T,N}}
        # allocation of `map` on tuple is optimized away
        checkoverflow && map(overflow_check, axes(parent), offsets)
        new{T, N, AA}(parent, offsets)
    end
end

"""
    OffsetVector(v, index)

Type alias and convenience constructor for one-dimensional [`OffsetArray`](@ref)s.
"""
const OffsetVector{T,AA<:AbstractVector{T}} = OffsetArray{T,1,AA}

"""
    OffsetMatrix(A, index1, index2)

Type alias and convenience constructor for two-dimensional [`OffsetArray`](@ref)s.
"""
const OffsetMatrix{T,AA<:AbstractMatrix{T}} = OffsetArray{T,2,AA}

# checks if the offset may be added to the range without overflowing
function overflow_check(r::AbstractUnitRange, offset::Integer)
    Base.hastypemax(eltype(r)) || return nothing
    # This gives some performance boost https://github.com/JuliaLang/julia/issues/33273
    throw_upper_overflow_error(val) = throw(OverflowError("offset should be <= $(typemax(Int) - val) corresponding to the axis $r, received an offset $offset"))
    throw_lower_overflow_error(val) = throw(OverflowError("offset should be >= $(typemin(Int) - val) corresponding to the axis $r, received an offset $offset"))

    # With ranges in the picture, first(r) might not necessarily be < last(r)
    # we therefore use the min and max of first(r) and last(r) to check for overflow
    firstlast_min, firstlast_max = minmax(first(r), last(r))

    if offset > 0 && firstlast_max > typemax(Int) - offset
        throw_upper_overflow_error(firstlast_max)
    elseif offset < 0 && firstlast_min < typemin(Int) - offset
        throw_lower_overflow_error(firstlast_min)
    end
    return nothing
end

# Tuples of integers are treated as offsets
# Empty Tuples are handled here
@inline function OffsetArray(A::AbstractArray, offsets::Tuple{Vararg{Integer}}; kw...)
    _checkindices(A, offsets, "offsets")
    OffsetArray{eltype(A), ndims(A), typeof(A)}(A, offsets; kw...)
end

# These methods are necessary to disallow incompatible dimensions for
# the OffsetVector and the OffsetMatrix constructors
for (FT, ND) in ((:OffsetVector, :1), (:OffsetMatrix, :2))
    @eval @inline function $FT(A::AbstractArray{<:Any,$ND}, offsets::Tuple{Vararg{Integer}}; kw...)
        _checkindices(A, offsets, "offsets")
        OffsetArray{eltype(A), $ND, typeof(A)}(A, offsets; kw...)
    end
    FTstr = string(FT)
    @eval @inline function $FT(A::AbstractArray, offsets::Tuple{Vararg{Integer}}; kw...)
        throw(ArgumentError($FTstr*" requires a "*string($ND)*"D array"))
    end
end

## OffsetArray constructors
for FT in (:OffsetArray, :OffsetVector, :OffsetMatrix)
    # Nested OffsetArrays may strip off the wrapper and collate the offsets
    # empty tuples are handled here
    @eval @inline function $FT(A::OffsetArray, offsets::Tuple{Vararg{Int}}; checkoverflow = true)
        _checkindices(A, offsets, "offsets")
        # ensure that the offsets may be added together without an overflow
        checkoverflow && map(overflow_check, axes(A), offsets)
        I = map(+, _offsets(A, parent(A)), offsets)
        $FT(parent(A), I, checkoverflow = false)
    end
    @eval @inline function $FT(A::OffsetArray, offsets::Tuple{Integer,Vararg{Integer}}; kw...)
        $FT(A, map(Int, offsets); kw...)
    end

    # In general, indices get converted to AbstractUnitRanges.
    # CartesianIndices{N} get converted to N ranges
    @eval @inline function $FT(A::AbstractArray, inds::Tuple{Any,Vararg{Any}}; kw...)
        $FT(A, _toAbstractUnitRanges(to_indices(A, axes(A), inds)); kw...)
    end

    # convert ranges to offsets
    @eval @inline function $FT(A::AbstractArray, inds::Tuple{AbstractUnitRange,Vararg{AbstractUnitRange}}; kw...)
        _checkindices(A, inds, "indices")
        # Performance gain by wrapping the error in a function: see https://github.com/JuliaLang/julia/issues/37558
        throw_dimerr(lA, lI) = throw(DimensionMismatch("supplied axes do not agree with the size of the array (got size $lA for the array and $lI for the indices"))
        lA = size(A)
        lI = map(length, inds)
        lA == lI || throw_dimerr(lA, lI)
        $FT(A, map(_offset, axes(A), inds); kw...)
    end

    @eval @inline $FT(A::AbstractArray, inds...; kw...) = $FT(A, inds; kw...)
    @eval @inline $FT(A::AbstractArray; checkoverflow = false) = $FT(A, ntuple(zero, Val(ndims(A))), checkoverflow = checkoverflow)

    @eval @inline $FT(A::AbstractArray, origin::Origin; checkoverflow = true) = $FT(A, origin(A); checkoverflow = checkoverflow)
end

# conversion-related methods
@inline OffsetArray{T}(M::AbstractArray, I...; kw...) where {T} = OffsetArray{T,ndims(M)}(M, I...; kw...)

@inline function OffsetArray{T,N}(M::AbstractArray{<:Any,N}, I...; kw...) where {T,N}
    M2 = _of_eltype(T, M)
    OffsetArray{T,N}(M2, I...; kw...)
end
@inline OffsetArray{T,N}(M::OffsetArray{T,N}, I...; kw...) where {T,N} = OffsetArray(M, I...; kw...)
@inline OffsetArray{T,N}(M::AbstractArray{T,N}, I...; kw...) where {T,N} = OffsetArray{T,N,typeof(M)}(M, I...; kw...)

@inline OffsetArray{T,N,A}(M::AbstractArray{<:Any,N}, I...; kw...) where {T,N,A<:AbstractArray{T,N}} = OffsetArray{T,N,A}(M, I; kw...)
@inline function OffsetArray{T,N,A}(M::AbstractArray{<:Any,N}, I::NTuple{N,Int}; checkoverflow = true) where {T,N,A<:AbstractArray{T,N}}
    checkoverflow && map(overflow_check, axes(M), I)
    Mv = no_offset_view(M)
    MvA = convert(A, Mv)::A
    Iof = map(+, _offsets(M), I)
    OffsetArray{T,N,A}(MvA, Iof, checkoverflow = false)
end
@inline function OffsetArray{T, N, AA}(parent::AbstractArray{<:Any,N}, offsets::NTuple{N, Integer}; kw...) where {T, N, AA<:AbstractArray{T,N}}
    OffsetArray{T, N, AA}(parent, map(Int, offsets)::NTuple{N,Int}; kw...)
end
@inline function OffsetArray{T,N,A}(M::AbstractArray{<:Any,N}, I::Tuple{AbstractUnitRange,Vararg{AbstractUnitRange}}; kw...) where {T,N,A<:AbstractArray{T,N}}
    _checkindices(M, I, "indices")
    # Performance gain by wrapping the error in a function: see https://github.com/JuliaLang/julia/issues/37558
    throw_dimerr(lA, lI) = throw(DimensionMismatch("supplied axes do not agree with the size of the array (got size $lA for the array and $lI for the indices"))
    lM = size(M)
    lI = map(length, I)
    lM == lI || throw_dimerr(lM, lI)
    OffsetArray{T,N,A}(M, map(_offset, axes(M), I); kw...)
end
@inline function OffsetArray{T,N,A}(M::AbstractArray{<:Any,N}, I::Tuple; kw...) where {T,N,A<:AbstractArray{T,N}}
    OffsetArray{T,N,A}(M, _toAbstractUnitRanges(to_indices(M, axes(M), I)); kw...)
end
@inline function OffsetArray{T,N,A}(M::AbstractArray{<:Any,N}; kw...) where {T,N,A<:AbstractArray{T,N}}
    Mv = no_offset_view(M)
    MvA = convert(A, Mv)::A
    OffsetArray{T,N,A}(MvA, _offsets(M); kw...)
end
@inline OffsetArray{T,N,A}(M::A; checkoverflow = false) where {T,N,A<:AbstractArray{T,N}} = OffsetArray{T,N,A}(M, ntuple(zero, Val(N)); checkoverflow = checkoverflow)

Base.convert(::Type{T}, M::AbstractArray) where {T<:OffsetArray} = M isa T ? M : T(M)

@inline AbstractArray{T,N}(M::OffsetArray{S,N}) where {T,S,N} = OffsetArray{T}(M)

# array initialization
@inline function OffsetArray{T,N}(init::ArrayInitializer, inds::Tuple{Vararg{OffsetAxisKnownLength}}; kw...) where {T,N}
    _checkindices(N, inds, "indices")
    AA = Array{T,N}(init, map(_indexlength, inds))
    OffsetArray{T, N, typeof(AA)}(AA, map(_indexoffset, inds); kw...)
end
@inline function OffsetArray{T, N}(init::ArrayInitializer, inds::Tuple; kw...) where {T, N}
    OffsetArray{T, N}(init, _toAbstractUnitRanges(inds); kw...)
end
@inline OffsetArray{T,N}(init::ArrayInitializer, inds...; kw...) where {T,N} = OffsetArray{T,N}(init, inds; kw...)

@inline OffsetArray{T}(init::ArrayInitializer, inds::NTuple{N, OffsetAxisKnownLength}; kw...) where {T,N} = OffsetArray{T,N}(init, inds; kw...)
@inline function OffsetArray{T}(init::ArrayInitializer, inds::Tuple; kw...) where {T}
    OffsetArray{T}(init, _toAbstractUnitRanges(inds); kw...)
end
@inline OffsetArray{T}(init::ArrayInitializer, inds...; kw...) where {T} = OffsetArray{T}(init, inds; kw...)

Base.IndexStyle(::Type{OA}) where {OA<:OffsetArray} = IndexStyle(parenttype(OA))
parenttype(::Type{OffsetArray{T,N,AA}}) where {T,N,AA} = AA
parenttype(A::OffsetArray) = parenttype(typeof(A))

Base.parent(A::OffsetArray) = A.parent

# TODO: Ideally we would delegate to the parent's broadcasting implementation, but that
#       is currently broken in sufficiently many implementation, namely RecursiveArrayTools, DistributedArrays
#       and StaticArrays, that it will take concentrated effort to get this working across the ecosystem.
#       The goal would be to have `OffsetArray(CuArray) .+ 1 == OffsetArray{CuArray}`.
# Base.Broadcast.BroadcastStyle(::Type{<:OffsetArray{<:Any, <:Any, AA}}) where AA = Base.Broadcast.BroadcastStyle(AA)

@inline Base.size(A::OffsetArray) = size(parent(A))

@inline Base.axes(A::OffsetArray) = map(IdOffsetRange, axes(parent(A)), A.offsets)
@inline Base.axes(A::OffsetArray, d) = d <= ndims(A) ? IdOffsetRange(axes(parent(A), d), A.offsets[d]) : IdOffsetRange(axes(parent(A), d))
@inline Base.axes1(A::OffsetArray{T,0}) where {T} = IdOffsetRange(axes(parent(A), 1))  # we only need to specialize this one

# Issue 128
# See https://github.com/JuliaLang/julia/issues/37274 for the issue reported in Base
# The fix https://github.com/JuliaLang/julia/pull/39404 should be available on v1.6
# The following method is added on older Julia versions to ensure correct behavior for OffsetVectors
if VERSION < v"1.6"
    @inline function Base.compute_linindex(A::OffsetVector, I::NTuple{N,Any}) where N
        IP = Base.fill_to_length(axes(A), Base.OneTo(1), Val(N))
        Base.compute_linindex(first(LinearIndices(A)), 1, IP, I)
    end
end

# Utils to translate a function to the parent while preserving offsets
unwrap(x) = x, identity
unwrap(x::OffsetArray) = parent(x), data -> OffsetArray(data, x.offsets, checkoverflow = false)
function parent_call(f, x)
    parent, wrap_offset = unwrap(x)
    wrap_offset(f(parent))
end

Base.similar(A::OffsetArray, ::Type{T}, dims::Dims) where T =
    similar(parent(A), T, dims)
function Base.similar(A::AbstractArray, ::Type{T}, shape::Tuple{OffsetAxisKnownLength,Vararg{OffsetAxisKnownLength}}) where T
    # strip IdOffsetRanges to extract the parent range and use it to generate the array
    new_shape = map(_strip_IdOffsetRange, shape)
    # route through _similar_axes_or_length to avoid a stack overflow if map(_strip_IdOffsetRange, shape) === shape
    # This tries to use new_shape directly in similar if similar(A, T, ::typeof(new_shape)) is defined
    # If this fails, it calls similar(A, T, map(_indexlength, new_shape)) to use the size along each axis
    # to generate the new array
    P = _similar_axes_or_length(A, T, new_shape, shape)
    return OffsetArray(P, map(_offset, axes(P), shape))
end
Base.similar(::Type{A}, sz::Tuple{Vararg{Int}}) where {A<:OffsetArray} = similar(Array{eltype(A)}, sz)
function Base.similar(::Type{T}, shape::Tuple{OffsetAxisKnownLength,Vararg{OffsetAxisKnownLength}}) where {T<:AbstractArray}
    new_shape = map(_strip_IdOffsetRange, shape)
    P = _similar_axes_or_length(T, new_shape, shape)
    OffsetArray(P, map(_offset, axes(P), shape))
end
# Try to use the axes to generate the parent array type
# This is useful if the axes have special meanings, such as with static arrays
# This method is hit if at least one axis provided to similar(A, T, axes) is an IdOffsetRange
# For example this is hit when similar(A::OffsetArray) is called,
# which expands to similar(A, eltype(A), axes(A))
_similar_axes_or_length(A, T, ax, ::Any) = similar(A, T, ax)
_similar_axes_or_length(AT, ax, ::Any) = similar(AT, ax)
# Handle the general case by resorting to lengths along each axis
# This is hit if none of the axes provided to similar(A, T, axes) are IdOffsetRanges,
# and if similar(A, T, axes::AX) is not defined for the type AX.
# In this case the best that we can do is to create a mutable array of the correct size
_similar_axes_or_length(A, T, ax::I, ::I) where {I} = similar(A, T, map(_indexlength, ax))
_similar_axes_or_length(AT, ax::I, ::I) where {I} = similar(AT, map(_indexlength, ax))

# reshape accepts a single colon
Base.reshape(A::AbstractArray, inds::OffsetAxis...) = reshape(A, inds)
function Base.reshape(A::AbstractArray, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}})
    AR = reshape(no_offset_view(A), map(_indexlength, inds))
    O = OffsetArray(AR, map(_offset, axes(AR), inds))
    return _popreshape(O, axes(AR), _filterreshapeinds(inds))
end

# Reshaping OffsetArrays can "pop" the original OffsetArray wrapper and return
# an OffsetArray(reshape(...)) instead of an OffsetArray(reshape(OffsetArray(...)))
# Short-circuit for AbstractVectors if the axes are compatible to get around the Base restriction
# to 1-based vectors
function _reshape(A::AbstractVector, inds::Tuple{OffsetAxis})
    @noinline throw_dimerr(ind::Integer) = throw(
        DimensionMismatch("parent has $(size(A,1)) elements, which is incompatible with length $ind"))
    @noinline throw_dimerr(ind) = throw(
        DimensionMismatch("parent has $(size(A,1)) elements, which is incompatible with indices $ind"))
    _checksize(first(inds), size(A,1)) || throw_dimerr(first(inds))
    A
end
_reshape(A, inds) = _reshape2(A, inds)
_reshape2(A, inds) = reshape(A, inds)
# avoid a stackoverflow by relegating to the parent if no_offset_view returns an offsetarray
_reshape2(A::OffsetArray, inds) = reshape(parent(A), inds)
_reshape_nov(A, inds) = _reshape(no_offset_view(A), inds)

Base.reshape(A::OffsetArray, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}}) =
    OffsetArray(_reshape(parent(A), inds), map(_toaxis, inds))
# And for non-offset axes, we can just return a reshape of the parent directly
Base.reshape(A::OffsetArray, inds::Tuple{Union{Integer,Base.OneTo},Vararg{Union{Integer,Base.OneTo}}}) = _reshape_nov(A, inds)
Base.reshape(A::OffsetArray, inds::Dims) = _reshape_nov(A, inds)
Base.reshape(A::OffsetVector, ::Colon) = A
Base.reshape(A::OffsetVector, ::Tuple{Colon}) = A
Base.reshape(A::OffsetArray, ::Colon) = reshape(A, (Colon(),))
Base.reshape(A::OffsetArray, inds::Union{Int,Colon}...) = reshape(A, inds)
Base.reshape(A::OffsetArray, inds::Tuple{Vararg{Union{Int,Colon}}}) = _reshape_nov(A, inds)

# permutedims in Base does not preserve axes, and can not be fixed in a non-breaking way
# This is a stopgap solution
Base.permutedims(v::OffsetVector) = reshape(v, (1, axes(v, 1)))

Base.fill(v, inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {N} =
    fill!(similar(Array{typeof(v)}, inds), v)
Base.zeros(::Type{T}, inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {T, N} =
    fill!(similar(Array{T}, inds), zero(T))
Base.ones(::Type{T}, inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {T, N} =
    fill!(similar(Array{T}, inds), one(T))
Base.trues(inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {N} =
    fill!(similar(BitArray, inds), true)
Base.falses(inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {N} =
    fill!(similar(BitArray, inds), false)

Base.zero(A::OffsetArray) = parent_call(zero, A)
Base.fill!(A::OffsetArray, x) = parent_call(Ap -> fill!(Ap, x), A)

## Indexing

# Note this gets the index of the parent *array*, not the index of the parent *range*
# Here's how one can think about this:
#   Δi = i - first(r)
#   i′ = first(r.parent) + Δi
# and one obtains the result below.
parentindex(r::IdOffsetRange, i) = i - r.offset

@propagate_inbounds Base.getindex(A::OffsetArray{<:Any,0})  = A.parent[]

@inline function Base.getindex(A::OffsetArray{<:Any,N}, I::Vararg{Int,N}) where N
    @boundscheck checkbounds(A, I...)
    J = map(parentindex, axes(A), I)
    @inbounds parent(A)[J...]
end

@propagate_inbounds Base.getindex(A::OffsetArray{<:Any,N}, c::Vararg{Colon,N}) where N =
    parent_call(x -> getindex(x, c...), A)

# With one Colon we use linear indexing.
# In this case we may forward the index to the parent, as the information about the axes is lost
# The exception to this is with OffsetVectors where the axis information is preserved,
# but that case is handled by getindex(::OffsetArray{<:Any,N}, ::Vararg{Colon,N})
@propagate_inbounds Base.getindex(A::OffsetArray, c::Colon) = A.parent[:]

@inline function Base.getindex(A::OffsetVector, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[parentindex(Base.axes1(A), i)]
end
@propagate_inbounds Base.getindex(A::OffsetArray, i::Int)  = parent(A)[i]

@inline function Base.setindex!(A::OffsetArray{T,N}, val, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    J = map(parentindex, axes(A), I)
    @inbounds parent(A)[J...] = val
    A
end

@inline function Base.setindex!(A::OffsetVector, val, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[parentindex(Base.axes1(A), i)] = val
    A
end
@propagate_inbounds function Base.setindex!(A::OffsetArray, val, i::Int)
    parent(A)[i] = val
    A
end

@inline Base.iterate(a::OffsetArray, i...) = iterate(parent(a), i...)

Base.in(x, A::OffsetArray) = in(x, parent(A))
Base.copy(A::OffsetArray) = parent_call(copy, A)

Base.strides(A::OffsetArray) = strides(parent(A))
Base.elsize(::Type{OffsetArray{T,N,A}}) where {T,N,A} = Base.elsize(A)
@inline Base.unsafe_convert(::Type{Ptr{T}}, A::OffsetArray{T}) where {T} = Base.unsafe_convert(Ptr{T}, parent(A))

# For fast broadcasting: ref https://discourse.julialang.org/t/why-is-there-a-performance-hit-on-broadcasting-with-offsetarrays/32194
Base.dataids(A::OffsetArray) = Base.dataids(parent(A))
Broadcast.broadcast_unalias(dest::OffsetArray, src::OffsetArray) = parent(dest) === parent(src) ? src : Broadcast.unalias(dest, src)

### Special handling for AbstractRange
const OffsetRange{T} = OffsetVector{T,<:AbstractRange{T}}
const OffsetUnitRange{T} = OffsetVector{T,<:AbstractUnitRange{T}}

Base.step(a::OffsetRange) = step(parent(a))

Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, or::OffsetRange) = Base.checkindex(Bool, inds, parent(or))

# Certain special methods for linear indexing with integer ranges (or OffsetRanges)
# These may bypass the default getindex(A, I...) pathway if the parent types permit this
# For example AbstractUnitRanges and Arrays have special linear indexing behavior defined

# If both the arguments are offset, we may unwrap the indices to call (::OffsetArray)[::AbstractRange{Int}]
@propagate_inbounds function Base.getindex(A::OffsetArray, r::OffsetRange{Int})
    _indexedby(A[parent(r)], axes(r))
end
# If the indices are offset, we may unwrap them and pass the parent to getindex
@propagate_inbounds function Base.getindex(A::AbstractRange, r::OffsetRange{Int})
    _indexedby(A[parent(r)], axes(r))
end

# An OffsetUnitRange might use the rapid getindex(::Array, ::AbstractUnitRange{Int}) for contiguous indexing
@propagate_inbounds function Base.getindex(A::Array, r::OffsetUnitRange{Int})
    B = A[_contiguousindexingtype(parent(r))]
    OffsetArray(B, axes(r), checkoverflow = false)
end

# avoid hitting the slow method getindex(::Array, ::AbstractRange{Int})
# instead use the faster getindex(::Array, ::UnitRange{Int})
if VERSION <= v"1.7.0-DEV.1039"
    @propagate_inbounds function Base.getindex(A::Array, r::Union{IdOffsetRange, IIUR})
        B = A[_contiguousindexingtype(r)]
        _indexedby(B, axes(r))
    end
end

# Linear Indexing of OffsetArrays with AbstractUnitRanges may use the faster contiguous indexing methods
@inline function Base.getindex(A::OffsetArray, r::AbstractUnitRange{Int})
    @boundscheck checkbounds(A, r)
    # nD OffsetArrays do not have their linear indices shifted, so we may forward the indices provided to the parent
    @inbounds B = parent(A)[_contiguousindexingtype(r)]
    _indexedby(B, axes(r))
end
@inline function Base.getindex(A::OffsetVector, r::AbstractUnitRange{Int})
    @boundscheck checkbounds(A, r)
    # OffsetVectors may have their linear indices shifted, so we subtract the offset from the indices provided
    @inbounds B = parent(A)[_subtractoffset(r, A.offsets[1])]
    _indexedby(B, axes(r))
end

# This method added mainly to index an OffsetRange with another range
@inline function Base.getindex(A::OffsetVector, r::AbstractRange{Int})
    @boundscheck checkbounds(A, r)
    @inbounds B = parent(A)[_subtractoffset(r, A.offsets[1])]
    _indexedby(B, axes(r))
end

# In general we would pass through getindex(A, I...) which calls to_indices(A, I) and finally to_index(I)
# An OffsetUnitRange{Int} has an equivalent IdOffsetRange with the same values and axes,
# something similar also holds for OffsetUnitRange{BigInt}
# We may replace the former with the latter in an indexing operation to obtain a performance boost
@inline function Base.to_index(r::OffsetUnitRange{<:Union{Int,BigInt}})
    of = first(axes(r,1)) - 1
    IdOffsetRange(_subtractoffset(parent(r), of), of)
end

@inline function _boundscheck_index_retaining_axes(r, s)
    @boundscheck checkbounds(r, s)
    @inbounds pr = r[UnitRange(s)]
    _indexedby(pr, axes(s))
end
@inline _boundscheck_return(r, s) = (@boundscheck checkbounds(r, s); s)

for OR in [:IIUR, :IdOffsetRange]
    for R in [:StepRange, :StepRangeLen, :LinRange, :UnitRange]
        @eval @inline Base.getindex(r::$R, s::$OR) = _boundscheck_index_retaining_axes(r, s)
    end

    # this method is needed for ambiguity resolution
    @eval @inline function Base.getindex(r::StepRangeLen{T,<:Base.TwicePrecision,<:Base.TwicePrecision}, s::$OR) where T
        _boundscheck_index_retaining_axes(r, s)
    end
end
Base.getindex(r::Base.OneTo, s::IdOffsetRange) = _boundscheck_index_retaining_axes(r, s)
if VERSION < v"1.7.0-beta2"
    Base.getindex(r::Base.OneTo, s::IIUR) = _boundscheck_index_retaining_axes(r, s)
end

# These methods are added to avoid ambiguities with Base.
# The ones involving Base types should be ported to Base and version-limited here
@inline Base.getindex(r::IdentityUnitRange, s::IIUR) = _boundscheck_return(r, s)
@inline Base.getindex(r::IdentityUnitRange, s::IdOffsetRange) = _boundscheck_return(r, s)
if IdentityUnitRange !== Base.Slice
    @inline Base.getindex(r::Base.Slice, s::IIUR) = _boundscheck_return(r, s)
    @inline Base.getindex(r::Base.Slice, s::IdOffsetRange) = _boundscheck_return(r, s)
end

# eltype conversion
# This may use specialized map methods for the parent
Base.map(::Type{T}, O::OffsetArray) where {T} = parent_call(x -> map(T, x), O)
Base.map(::Type{T}, r::IdOffsetRange) where {T<:Real} = _indexedby(map(T, UnitRange(r)), axes(r))
if eltype(IIUR) === Int
    # This is type-piracy, but there is no way to convert an IdentityUnitRange to a non-Int type in Base
    Base.map(::Type{T}, r::IdentityUnitRange) where {T<:Real} = _indexedby(map(T, UnitRange(r)), axes(r))
end

# mapreduce is faster with an IdOffsetRange than with an OffsetUnitRange
# We therefore convert OffsetUnitRanges to IdOffsetRanges with the same values and axes
function Base.mapreduce(f, op, A1::OffsetUnitRange{<:Integer}, As::OffsetUnitRange{<:Integer}...; kw...)
    As = (A1, As...)
    ofs = map(A -> first(axes(A,1)) - 1, As)
    AIds = map((A, of) -> IdOffsetRange(_subtractoffset(parent(A), of), of), As, ofs)
    mapreduce(f, op, AIds...; kw...)
end

# Optimize certain reductions that treat an OffsetVector as a list
for f in [:minimum, :maximum, :extrema, :sum]
    @eval Base.$f(r::OffsetRange) = $f(parent(r))
end

function Base.show(io::IO, r::OffsetRange)
    show(io, r.parent)
    print(io, " with indices ", UnitRange(axes(r, 1)))
end
Base.show(io::IO, ::MIME"text/plain", r::OffsetRange) = show(io, r)


### Some mutating functions defined only for OffsetVector ###

Base.resize!(A::OffsetVector, nl::Integer) = (resize!(A.parent, nl); A)
Base.push!(A::OffsetVector, x...) = (push!(A.parent, x...); A)
Base.pop!(A::OffsetVector) = pop!(A.parent)
Base.append!(A::OffsetVector, items) = (append!(A.parent, items); A)
Base.empty!(A::OffsetVector) = (empty!(A.parent); A)

# These functions keep the summary compact
function Base.inds2string(inds::Tuple{Vararg{Union{IdOffsetRange, IdentityUnitRange{<:IdOffsetRange}}}})
    Base.inds2string(map(UnitRange, inds))
end
Base.showindices(io::IO, ind1::IdOffsetRange, inds::IdOffsetRange...) = Base.showindices(io, map(UnitRange, (ind1, inds...))...)

function Base.showarg(io::IO, @nospecialize(a::OffsetArray), toplevel)
    print(io, "OffsetArray(")
    Base.showarg(io, parent(a), false)
    Base.showindices(io, axes(a)...)
    print(io, ')')
    if toplevel
        print(io, " with eltype ", eltype(a))
    end
end

function Base.replace_in_print_matrix(A::OffsetArray{<:Any,2}, i::Integer, j::Integer, s::AbstractString)
    J = map(parentindex, axes(A), (i,j))
    Base.replace_in_print_matrix(parent(A), J..., s)
end
function Base.replace_in_print_matrix(A::OffsetArray{<:Any,1}, i::Integer, j::Integer, s::AbstractString)
    ip = parentindex(axes(A,1), i)
    Base.replace_in_print_matrix(parent(A), ip, j, s)
end

"""
    no_offset_view(A)

Return an `AbstractArray` that shares structure and underlying data with the argument,
but uses 1-based indexing. May just return the argument when applicable.
Not exported.

The default implementation uses `OffsetArrays`, but other types should use something more
specific to remove a level of indirection when applicable.

```jldoctest; setup=:(using OffsetArrays)
julia> A = [1 3 5; 2 4 6];

julia> O = OffsetArray(A, 0:1, -1:1)
2×3 OffsetArray(::$(Matrix{Int}), 0:1, -1:1) with eltype $Int with indices 0:1×-1:1:
 1  3  5
 2  4  6

julia> OffsetArrays.no_offset_view(O)[1,1] = -9
-9

julia> A
2×3 $(Matrix{Int}):
 -9  3  5
  2  4  6
```
"""
no_offset_view(A::OffsetArray) = no_offset_view(parent(A))
if isdefined(Base, :IdentityUnitRange)
    # valid only if Slice is distinguished from IdentityUnitRange
    no_offset_view(a::Base.Slice{<:Base.OneTo}) = a
    no_offset_view(a::Base.Slice) = Base.Slice(UnitRange(a))
    no_offset_view(S::SubArray) = view(parent(S), map(no_offset_view, parentindices(S))...)
end
no_offset_view(a::Array) = a
no_offset_view(i::Number) = i
no_offset_view(A::AbstractArray) = _no_offset_view(axes(A), A)
_no_offset_view(::Tuple{}, A::AbstractArray{T,0}) where T = A
_no_offset_view(::Tuple{Base.OneTo, Vararg{Base.OneTo}}, A::AbstractArray) = A
# the following method is needed for ambiguity resolution
_no_offset_view(::Tuple{Base.OneTo, Vararg{Base.OneTo}}, A::AbstractUnitRange) = A
_no_offset_view(::Any, A::AbstractArray) = OffsetArray(A, Origin(1))
_no_offset_view(::Any, A::AbstractUnitRange) = UnitRange(A)

#####
# center/centered
# These two helpers are deliberately not exported; their meaning can be very different in
# other scenarios and will be very likely to cause name conflicts if exported.
#####
"""
    center(A, [r::RoundingMode=RoundDown])::Dims

Return the center coordinate of given array `A`. If `size(A, k)` is even,
a rounding procedure will be applied with mode `r`.

!!! compat "OffsetArrays 1.9"
    This method requires at least OffsetArrays 1.9.

# Examples

```jldoctest; setup=:(using OffsetArrays)
julia> A = reshape(collect(1:9), 3, 3)
3×3 $(Matrix{Int}):
 1  4  7
 2  5  8
 3  6  9

julia> c = OffsetArrays.center(A)
(2, 2)

julia> A[c...]
5

julia> Ao = OffsetArray(A, -2, -2); # axes (-1:1, -1:1)

julia> c = OffsetArrays.center(Ao)
(0, 0)

julia> Ao[c...]
5
```

To shift the center coordinate of the given array to `(0, 0, ...)`, you
can use [`centered`](@ref OffsetArrays.centered).
"""
function center(A::AbstractArray, r::RoundingMode=RoundDown)
    map(axes(A)) do inds
        round(Int, (length(inds)-1)/2, r) + first(inds)
    end
end

"""
    centered(A, cp=center(A)) -> Ao

Shift the center coordinate/point `cp` of array `A` to `(0, 0, ..., 0)`. Internally, this is
equivalent to `OffsetArray(A, .-cp)`.

!!! compat "OffsetArrays 1.9"
    This method requires at least OffsetArrays 1.9.

# Examples

```jldoctest; setup=:(using OffsetArrays)
julia> A = reshape(collect(1:9), 3, 3)
3×3 $(Matrix{Int}):
 1  4  7
 2  5  8
 3  6  9

julia> Ao = OffsetArrays.centered(A); # axes (-1:1, -1:1)

julia> Ao[0, 0]
5

julia> Ao = OffsetArray(A, OffsetArrays.Origin(0)); # axes (0:2, 0:2)

julia> Aoo = OffsetArrays.centered(Ao); # axes (-1:1, -1:1)

julia> Aoo[0, 0]
5
```

Users are allowed to pass `cp` to change how "center point" is interpreted, but the meaning of the
output array should be reinterpreted as well. For instance, if `cp = map(last, axes(A))` then this
function no longer shifts the center point but instead the bottom-right point to `(0, 0, ..., 0)`.
A commonly usage of `cp` is to change the rounding behavior when the array is of even size at some
dimension:

```jldoctest; setup=:(using OffsetArrays)
julia> A = reshape(collect(1:4), 2, 2) # Ideally the center should be (1.5, 1.5) but OffsetArrays only support integer offsets
2×2 $(Matrix{Int}):
 1  3
 2  4

julia> OffsetArrays.centered(A, OffsetArrays.center(A, RoundUp)) # set (2, 2) as the center point
2×2 OffsetArray(::$(Matrix{Int}), -1:0, -1:0) with eltype $(Int) with indices -1:0×-1:0:
 1  3
 2  4

julia> OffsetArrays.centered(A, OffsetArrays.center(A, RoundDown)) # set (1, 1) as the center point
2×2 OffsetArray(::$(Matrix{Int}), 0:1, 0:1) with eltype $(Int) with indices 0:1×0:1:
 1  3
 2  4
```

See also [`center`](@ref OffsetArrays.center).
"""
centered(A::AbstractArray, cp::Dims=center(A)) = OffsetArray(A, .-cp)

centered(A::AbstractArray, i::CartesianIndex) = centered(A, Tuple(i))

####
# work around for segfault in searchsorted*
#  https://github.com/JuliaLang/julia/issues/33977
####

function _safe_searchsorted(v::OffsetArray, x, ilo::T, ihi::T, o::Base.Ordering) where T<:Integer
    u = T(1)
    lo = ilo - u
    hi = ihi + u
    @inbounds while lo < hi - u
        m = (lo + hi) ÷ 2
        if Base.lt(o, v[m], x)
            lo = m
        elseif Base.lt(o, x, v[m])
            hi = m
        else
            a = searchsortedfirst(v, x, max(lo,ilo), m, o)
            b = searchsortedlast(v, x, m, min(hi,ihi), o)
            return a : b
        end
    end
    return (lo + 1) : (hi - 1)
end
function _safe_searchsortedfirst(v::OffsetArray, x, lo::T, hi::T, o::Base.Ordering) where T<:Integer
    u = T(1)
    lo = lo - u
    hi = hi + u
    @inbounds while lo < hi - u
        m = (lo + hi) ÷ 2
        if Base.lt(o, v[m], x)
            lo = m
        else
            hi = m
        end
    end
    return hi
end
function _safe_searchsortedlast(v::OffsetArray, x, lo::T, hi::T, o::Base.Ordering) where T<:Integer
    u = T(1)
    lo = lo - u
    hi = hi + u
    @inbounds while lo < hi - u
        m = (lo + hi) ÷ 2
        if Base.lt(o, x, v[m])
            hi = m
        else
            lo = m
        end
    end
    return lo
end

if VERSION ≤ v"1.2"
    # ambiguity warnings in earlier versions
    Base.searchsorted(v::OffsetArray, x, ilo::Int, ihi::Int, o::Base.Ordering) =
        _safe_searchsorted(v, x, ilo, ihi, o)
    Base.searchsortedfirst(v::OffsetArray, x, lo::Int, hi::Int, o::Base.Ordering) =
        _safe_searchsortedfirst(v, x, lo, hi, o)
    Base.searchsortedlast(v::OffsetArray, x, lo::Int, hi::Int, o::Base.Ordering) =
        _safe_searchsortedlast(v, x, lo, hi, o)
end

Base.searchsorted(v::OffsetArray, x, ilo::T, ihi::T, o::Base.Ordering) where T<:Integer =
    _safe_searchsorted(v, x, ilo, ihi, o)
Base.searchsortedfirst(v::OffsetArray, x, lo::T, hi::T, o::Base.Ordering) where T<:Integer =
    _safe_searchsortedfirst(v, x, lo, hi, o)
Base.searchsortedlast(v::OffsetArray, x, lo::T, hi::T, o::Base.Ordering) where T<:Integer =
    _safe_searchsortedlast(v, x, lo, hi, o)

if VERSION < v"1.1.0-DEV.783"
    Base.copyfirst!(dest::OffsetArray, src::OffsetArray) = (maximum!(parent(dest), parent(src)); return dest)
end

if VERSION <= v"1.7.0-DEV.400"
    # https://github.com/JuliaLang/julia/pull/39393
    # index for zero-argument getindex should be first linear index instead of 1 (#194)
    Base._to_linear_index(A::OffsetArray) = first(LinearIndices(A))
end

##
# Adapt allows for automatic conversion of CPU OffsetArrays to GPU OffsetArrays
##
import Adapt
Adapt.adapt_structure(to, O::OffsetArray) = parent_call(x -> Adapt.adapt(to, x), O)

if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
    _precompile_()
end


##
# Deprecations
##

# This is a bad API design as it introduces counter intuitive results (#250)
@deprecate centered(A::AbstractArray, r::RoundingMode) OffsetArray(A, .-center(A, r)) false

end # module
