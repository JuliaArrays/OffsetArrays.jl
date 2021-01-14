module OffsetArrays

using Base: tail, @propagate_inbounds
@static if !isdefined(Base, :IdentityUnitRange)
    const IdentityUnitRange = Base.Slice
else
    using Base: IdentityUnitRange
end

export OffsetArray, OffsetMatrix, OffsetVector

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
    function OffsetArray{T, N, AA}(parent::AA, offsets::NTuple{N, Int}) where {T, N, AA<:AbstractArray{T,N}}
        # allocation of `map` on tuple is optimized away
        map(overflow_check, axes(parent), offsets)
        new{T, N, AA}(parent, offsets)
    end
end

function OffsetArray{T, N, AA}(parent::AA, offsets::NTuple{N, <:Integer}) where {T, N, AA<:AbstractArray{T,N}}
    OffsetArray{T, N, AA}(parent, map(x -> convert(Int, x), offsets))
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

function overflow_check(r, offset::T) where T
    # This gives some performance boost https://github.com/JuliaLang/julia/issues/33273
    throw_upper_overflow_error() = throw(OverflowError("Boundary overflow detected: offset should be <= $(typemax(T) - last(r)) for offsets of type $T, received $offset"))
    throw_lower_overflow_error() = throw(OverflowError("Boundary overflow detected: offset should be >= $(typemin(T) - first(r)) for offsets of type $T, received $offset"))

    if offset > 0 && last(r) > typemax(T) - offset
        throw_upper_overflow_error()
    elseif offset < 0 && first(r) < typemin(T) - offset
        throw_lower_overflow_error()
    end
end

# Tuples of integers are treated as offsets
# Empty Tuples are handled here
@inline function OffsetArray(A::AbstractArray, offsets::Tuple{Vararg{Integer}})
    _checkindices(A, offsets, "offsets")
    OffsetArray{eltype(A), ndims(A), typeof(A)}(A, offsets)
end

# These methods are necessary to disallow incompatible dimensions for
# the OffsetVector and the OffsetMatrix constructors
for (FT, ND) in ((:OffsetVector, :1), (:OffsetMatrix, :2))
    @eval @inline function $FT(A::AbstractArray{<:Any,$ND}, offsets::Tuple{Vararg{Integer}})
        _checkindices(A, offsets, "offsets")
        OffsetArray{eltype(A), $ND, typeof(A)}(A, offsets)
    end
    FTstr = string(FT)
    @eval @inline function $FT(A::AbstractArray, offsets::Tuple{Vararg{Integer}})
        throw(ArgumentError($FTstr*" requires a "*string($ND)*"D array"))
    end
end

## OffsetArray constructors
for FT in (:OffsetArray, :OffsetVector, :OffsetMatrix)
    # Nested OffsetArrays may strip off the wrapper and collate the offsets
    @eval @inline function $FT(A::OffsetArray, offsets::Tuple{Vararg{Integer}})
        _checkindices(A, offsets, "offsets")
        $FT(parent(A), map(+, A.offsets, offsets))
    end

    # In general, indices get converted to AbstractUnitRanges.
    # CartesianIndices{N} get converted to N ranges
    @eval @inline function $FT(A::AbstractArray, inds::Tuple{Any,Vararg{Any}})
        $FT(A, _toAbstractUnitRanges(to_indices(A, axes(A), inds)))
    end

    # convert ranges to offsets
    @eval @inline function $FT(A::AbstractArray, inds::Tuple{AbstractUnitRange,Vararg{AbstractUnitRange}})
        _checkindices(A, inds, "indices")
        # Performance gain by wrapping the error in a function: see https://github.com/JuliaLang/julia/issues/37558
        throw_dimerr(lA, lI) = throw(DimensionMismatch("supplied axes do not agree with the size of the array (got size $lA for the array and $lI for the indices"))
        lA = size(A)
        lI = map(length, inds)
        lA == lI || throw_dimerr(lA, lI)
        $FT(A, map(_offset, axes(A), inds))
    end

    @eval @inline $FT(A::AbstractArray, inds::Vararg) = $FT(A, inds)

    @eval @inline $FT(A::AbstractArray, origin::Origin) = $FT(A, origin(A))
end

# array initialization
@inline function OffsetArray{T,N}(init::ArrayInitializer, inds::Tuple{Vararg{OffsetAxisKnownLength}}) where {T,N}
    _checkindices(N, inds, "indices")
    AA = Array{T,N}(init, map(_indexlength, inds))
    OffsetArray{T, N, typeof(AA)}(AA, map(_indexoffset, inds))
end
@inline function OffsetArray{T, N}(init::ArrayInitializer, inds::Tuple) where {T, N}
    OffsetArray{T, N}(init, _toAbstractUnitRanges(inds))
end
@inline OffsetArray{T,N}(init::ArrayInitializer, inds::Vararg) where {T,N} = OffsetArray{T,N}(init, inds)

@inline OffsetArray{T}(init::ArrayInitializer, inds::NTuple{N, OffsetAxisKnownLength}) where {T,N} = OffsetArray{T,N}(init, inds)
@inline function OffsetArray{T}(init::ArrayInitializer, inds::Tuple) where {T}
    OffsetArray{T}(init, _toAbstractUnitRanges(inds))
end
@inline OffsetArray{T}(init::ArrayInitializer, inds::Vararg) where {T} = OffsetArray{T}(init, inds)

Base.IndexStyle(::Type{OA}) where {OA<:OffsetArray} = IndexStyle(parenttype(OA))
parenttype(::Type{OffsetArray{T,N,AA}}) where {T,N,AA} = AA
parenttype(A::OffsetArray) = parenttype(typeof(A))

Base.parent(A::OffsetArray) = A.parent

# TODO: Ideally we would delegate to the parent's broadcasting implementation, but that
#       is currently broken in sufficiently many implementation, namely RecursiveArrayTools, DistributedArrays
#       and StaticArrays, that it will take concentrated effort to get this working across the ecosystem.
#       The goal would be to have `OffsetArray(CuArray) .+ 1 == OffsetArray{CuArray}`.
# Base.Broadcast.BroadcastStyle(::Type{<:OffsetArray{<:Any, <:Any, AA}}) where AA = Base.Broadcast.BroadcastStyle(AA)

Base.eachindex(::IndexCartesian, A::OffsetArray) = CartesianIndices(axes(A))
Base.eachindex(::IndexLinear, A::OffsetVector)   = axes(A, 1)

@inline Base.size(A::OffsetArray) = size(parent(A))
@inline Base.size(A::OffsetArray, d) = size(parent(A), d)

@inline Base.axes(A::OffsetArray) = map(IdOffsetRange, axes(parent(A)), A.offsets)
@inline Base.axes(A::OffsetArray, d) = d <= ndims(A) ? IdOffsetRange(axes(parent(A), d), A.offsets[d]) : IdOffsetRange(axes(parent(A), d))
@inline Base.axes1(A::OffsetArray{T,0}) where {T} = IdOffsetRange(axes(parent(A), 1))  # we only need to specialize this one

Base.similar(A::OffsetArray, ::Type{T}, dims::Dims) where T =
    similar(parent(A), T, dims)
function Base.similar(A::AbstractArray, ::Type{T}, inds::Tuple{OffsetAxisKnownLength,Vararg{OffsetAxisKnownLength}}) where T
    B = similar(A, T, map(_indexlength, inds))
    return OffsetArray(B, map(_offset, axes(B), inds))
end

# reshape accepts a single colon
Base.reshape(A::AbstractArray, inds::OffsetAxis...) = reshape(A, inds)
function Base.reshape(A::AbstractArray, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}})
    AR = reshape(A, map(_indexlength, inds))
    return OffsetArray(AR, map(_offset, axes(AR), inds))
end

# Reshaping OffsetArrays can "pop" the original OffsetArray wrapper and return
# an OffsetArray(reshape(...)) instead of an OffsetArray(reshape(OffsetArray(...)))
Base.reshape(A::OffsetArray, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}}) =
    OffsetArray(reshape(parent(A), map(_indexlength, inds)), map(_indexoffset, inds))
# And for non-offset axes, we can just return a reshape of the parent directly
Base.reshape(A::OffsetArray, inds::Tuple{Union{Integer,Base.OneTo},Vararg{Union{Integer,Base.OneTo}}}) = reshape(parent(A), inds)
Base.reshape(A::OffsetArray, inds::Dims) = reshape(parent(A), inds)
Base.reshape(A::OffsetArray, ::Colon) = reshape(parent(A), Colon())
Base.reshape(A::OffsetVector, ::Colon) = A
Base.reshape(A::OffsetArray, inds::Union{Int,Colon}...) = reshape(parent(A), inds)
Base.reshape(A::OffsetArray, inds::Tuple{Vararg{Union{Int,Colon}}}) = reshape(parent(A), inds)

function Base.similar(::Type{T}, shape::Tuple{OffsetAxisKnownLength,Vararg{OffsetAxisKnownLength}}) where {T<:AbstractArray}
    P = T(undef, map(_indexlength, shape))
    OffsetArray(P, map(_offset, axes(P), shape))
end

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

## Indexing

# Note this gets the index of the parent *array*, not the index of the parent *range*
# Here's how one can think about this:
#   Δi = i - first(r)
#   i′ = first(r.parent) + Δi
# and one obtains the result below.
parentindex(r::IdOffsetRange, i) = i - r.offset

@inline function Base.getindex(A::OffsetArray{T,N}, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    J = map(parentindex, axes(A), I)
    @inbounds parent(A)[J...]
end

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

Base.strides(A::OffsetArray) = strides(parent(A))
Base.elsize(::Type{OffsetArray{T,N,A}}) where {T,N,A} = Base.elsize(A)
@inline Base.unsafe_convert(::Type{Ptr{T}}, A::OffsetArray{T}) where {T} = Base.unsafe_convert(Ptr{T}, parent(A))

# For fast broadcasting: ref https://discourse.julialang.org/t/why-is-there-a-performance-hit-on-broadcasting-with-offsetarrays/32194
Base.dataids(A::OffsetArray) = Base.dataids(parent(A))
Broadcast.broadcast_unalias(dest::OffsetArray, src::OffsetArray) = parent(dest) === parent(src) ? src : Broadcast.unalias(dest, src)

### Special handling for AbstractRange

const OffsetRange{T} = OffsetArray{T,1,<:AbstractRange{T}}
const IIUR = IdentityUnitRange{S} where S<:AbstractUnitRange{T} where T<:Integer

Base.step(a::OffsetRange) = step(parent(a))

@propagate_inbounds Base.getindex(a::OffsetRange, r::OffsetRange) = OffsetArray(a[parent(r)], r.offsets)
@propagate_inbounds function Base.getindex(a::OffsetRange, r::IdOffsetRange)
    OffsetArray(a.parent[r.parent .+ (r.offset - a.offsets[1])], r.offset)
end
@propagate_inbounds Base.getindex(r::OffsetRange, s::IIUR) =
    OffsetArray(r[s.indices], s)
@propagate_inbounds Base.getindex(a::OffsetRange, r::AbstractRange) = a.parent[r .- a.offsets[1]]
@propagate_inbounds Base.getindex(a::AbstractRange, r::OffsetRange) = OffsetArray(a[parent(r)], r.offsets)

for OR in [:IIUR, :IdOffsetRange]
    for R in [:StepRange, :StepRangeLen, :LinRange, :UnitRange]
        @eval @propagate_inbounds Base.getindex(r::$R, s::$OR) = OffsetArray(r[_unwrap(s)], axes(s))
    end

    # this method is needed for ambiguity resolution
    @eval @propagate_inbounds Base.getindex(r::StepRangeLen{T,<:Base.TwicePrecision,<:Base.TwicePrecision}, s::$OR) where T =
    OffsetArray(r[_unwrap(s)], axes(s))
end

#= Integer UnitRanges may return an appropriate AbstractUnitRange{<:Integer}, as the result may be used in indexing, and
indexing is faster with ranges =#
@propagate_inbounds Base.getindex(r::UnitRange{<:Integer}, s::IdOffsetRange) = IdOffsetRange(r[_unwrap(s)] .- s.offset, s.offset)
@propagate_inbounds Base.getindex(r::UnitRange{<:Integer}, s::IIUR) = IdentityUnitRange(r[_unwrap(s)])

function Base.show(io::IO, r::OffsetRange)
    show(io, r.parent)
    o = r.offsets[1]
    print(io, " with indices ", o+1:o+length(r))
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
function no_offset_view(A::AbstractArray)
    if Base.has_offset_axes(A)
        OffsetArray(A, Origin(1))
    else
        A
    end
end

no_offset_view(A::OffsetArray) = no_offset_view(parent(A))
no_offset_view(a::Array) = a

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

##
# Adapt allows for automatic conversion of CPU OffsetArrays to GPU OffsetArrays
##
import Adapt
Adapt.adapt_structure(to, x::OffsetArray) = OffsetArray(Adapt.adapt(to, parent(x)), x.offsets)

if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
    _precompile_()
end

end # module
