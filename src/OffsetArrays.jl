module OffsetArrays

using Base: Indices, tail, @propagate_inbounds
@static if !isdefined(Base, :IdentityUnitRange)
    const IdentityUnitRange = Base.Slice
else
    using Base: IdentityUnitRange
end

export OffsetArray, OffsetMatrix, OffsetVector

include("axes.jl")
include("utils.jl")

# Technically we know the length of CartesianIndices but we need to convert it first, so here we
# don't put it in OffsetAxisKnownLength.
const OffsetAxisKnownLength = Union{Integer, AbstractUnitRange, IdOffsetRange}
const OffsetAxis = Union{OffsetAxisKnownLength, CartesianIndices, Colon}
const ArrayInitializer = Union{UndefInitializer, Missing, Nothing}

## OffsetArray
"""
    OffsetArray(A, indices...)

Return an `AbstractArray` that shares element type and size with the first argument, but
used the given `indices`, which are checked for compatible size.

# Example

There are two types of `indices`: integers and ranges-like types.

Integers are recognized as offsets, where `0` means no offsets are applied:

```jldoctest; setup=:(using OffsetArrays)
julia> A = OffsetArray(reshape(1:6, 2, 3), -1, -2)
2×3 OffsetArray(reshape(::UnitRange{Int64}, 2, 3), 0:1, -1:1) with eltype Int64 with indices 0:1×-1:1:
 1  3  5
 2  4  6

julia> A[0, 1]
5
```

Examples of range-like types are: `Colon()`(aka `:`), `UnitRange`(e.g, `-1:2`), and
`CartesianIndices`.

```jldoctest; setup=:(using OffsetArrays)
julia> OffsetArray(reshape(1:6, 2, 3), 0:1, -1:1)
2×3 OffsetArray(reshape(::UnitRange{Int64}, 2, 3), 0:1, -1:1) with eltype Int64 with indices 0:1×-1:1:
 1  3  5
 2  4  6

julia> OffsetArray(reshape(1:6, 2, 3), :, -1:1) # : as a placeholder means no offset is applied at this dimension
2×3 OffsetArray(reshape(::UnitRange{Int64}, 2, 3), 1:2, -1:1) with eltype Int64 with indices 1:2×-1:1:
 1  3  5
 2  4  6

julia> OffsetArray(reshape(1:6, 2, 3), CartesianIndex(0, -1):CartesianIndex(1, 1))
2×3 OffsetArray(reshape(::UnitRange{Int64}, 2, 3), 0:1, -1:1) with eltype Int64 with indices 0:1×-1:1:
 1  3  5
 2  4  6
```

Integers and range-like types can't be used interchangebly:

```julia
julia> OffsetArray(reshape(1:6, 2, 3), 0, -1:1)
ERROR: [...]
```

"""
struct OffsetArray{T,N,AA<:AbstractArray} <: AbstractArray{T,N}
    parent::AA
    offsets::NTuple{N,Int}
    function OffsetArray{T, N, AA}(parent::AA, offsets::NTuple{N, Int}) where {T, N, AA<:AbstractArray}
        @boundscheck overflow_check.(axes(parent), offsets)
        new{T, N, AA}(parent, offsets)
    end
end
const OffsetVector{T,AA<:AbstractArray} = OffsetArray{T,1,AA}
const OffsetMatrix{T,AA<:AbstractArray} = OffsetArray{T,2,AA}

function overflow_check(r, offset::T) where T
    # This gives some performance boost https://github.com/JuliaLang/julia/issues/33273
    throw_upper_overflow_error() = throw(ArgumentError("Boundary overflow detected: offset $offset should be equal or less than $(typemax(T) - last(r))"))
    throw_lower_overflow_error() = throw(ArgumentError("Boundary overflow detected: offset $offset should be equal or greater than $(typemin(T) - first(r))"))

    if offset > 0 && last(r) > typemax(T) - offset
        throw_upper_overflow_error()
    elseif offset < 0 && first(r) < typemin(T) - offset
        throw_lower_overflow_error()
    end
end
## OffsetArray constructors

for FT in (:OffsetArray, :OffsetVector, :OffsetMatrix)
    # The only route out to inner constructor
    @eval function $FT(A::AbstractArray{T, N}, offsets::NTuple{N, Integer}) where {T, N}
        ndims(A) == N || throw(DimensionMismatch("The number of offsets $(N) should equal ndims(A) = $(ndims(A))"))
        OffsetArray{T, ndims(A), typeof(A)}(A, offsets)
    end
    # nested OffsetArrays
    @eval $FT(A::OffsetArray{T, N}, offsets::NTuple{N, Integer}) where {T,N} = $FT(parent(A), A.offsets .+ offsets)
    # convert ranges to offsets
    @eval function $FT(A::AbstractArray{T}, inds::NTuple{N,OffsetAxisKnownLength}) where {T,N}
        axparent = axes(A)
        lA = map(length, axparent)
        lI = map(length, inds)
        lA == lI || throw(DimensionMismatch("supplied axes do not agree with the size of the array (got size $lA for the array and $lI for the indices"))
        $FT(A, map(_offset, axparent, inds)) 
    end
    # lower CartesianIndices and Colon
    @eval function $FT(A::AbstractArray{T}, inds::NTuple{N, OffsetAxis}) where {T, N}
        indsN = _uncolonindices(A, _expandCartesianIndices(inds))
        $FT(A, indsN)
    end
    @eval $FT(A::AbstractArray{T}, inds::Vararg{OffsetAxis,N}) where {T, N} = $FT(A, inds)
end

# array initialization
function OffsetArray{T,N}(init::ArrayInitializer, inds::NTuple{N, OffsetAxisKnownLength}) where {T,N}
    AA = Array{T,N}(init, map(_indexlength, inds))
    OffsetArray{T, N, typeof(AA)}(AA, map(_indexoffset, inds))
end
function OffsetArray{T, N}(init::ArrayInitializer, inds::NTuple{NT, Union{OffsetAxisKnownLength, CartesianIndices}}) where {T, N, NT}
    # NT is probably not the actual dimension of the array; CartesianIndices might contain multiple dimensions
    indsN = _expandCartesianIndices(inds)
    length(indsN) == N || throw(DimensionMismatch("The number of offsets $(length(indsN)) should equal ndims(A) = $N"))
    OffsetArray{T, N}(init, indsN)
end
OffsetArray{T,N}(init::ArrayInitializer, inds::Union{OffsetAxisKnownLength, CartesianIndices}...) where {T,N} = OffsetArray{T,N}(init, inds)

OffsetArray{T}(init::ArrayInitializer, inds::NTuple{N, OffsetAxisKnownLength}) where {T,N} = OffsetArray{T,N}(init, inds)
function OffsetArray{T}(init::ArrayInitializer, inds::NTuple{N, Union{OffsetAxisKnownLength, CartesianIndices}}) where {T, N}
    # N is probably not the actual dimension of the array; CartesianIndices might contain multiple dimensions
    indsN = _expandCartesianIndices(inds)
    OffsetArray{T, length(indsN)}(init, indsN)
end
OffsetArray{T}(init::ArrayInitializer, inds::Union{OffsetAxisKnownLength, CartesianIndices}...) where {T} = OffsetArray{T}(init, inds)

Base.IndexStyle(::Type{OA}) where {OA<:OffsetArray} = IndexStyle(parenttype(OA))
parenttype(::Type{OffsetArray{T,N,AA}}) where {T,N,AA} = AA
parenttype(A::OffsetArray) = parenttype(typeof(A))

Base.parent(A::OffsetArray) = A.parent

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

function Base.similar(::Type{T}, shape::Tuple{OffsetAxis,Vararg{OffsetAxis}}) where {T<:AbstractArray}
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
    J = @inbounds map(parentindex, axes(A), I)
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

@propagate_inbounds Base.getindex(r::UnitRange, s::IIUR) =
    OffsetArray(r[s.indices], s)

@propagate_inbounds Base.getindex(r::StepRange, s::IIUR) =
    OffsetArray(r[s.indices], s)

# this method is needed for ambiguity resolution
@propagate_inbounds Base.getindex(r::StepRangeLen{T,<:Base.TwicePrecision,<:Base.TwicePrecision}, s::IIUR) where T =
    OffsetArray(r[s.indices], s)

@propagate_inbounds Base.getindex(r::StepRangeLen{T}, s::IIUR) where {T} =
    OffsetArray(r[s.indices], s)

@propagate_inbounds Base.getindex(r::LinRange, s::IIUR) =
    OffsetArray(r[s.indices], s)

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
    
function Base.showarg(io::IO, a::OffsetArray, toplevel)
    print(io, "OffsetArray(")
    Base.showarg(io, parent(a), false)
    Base.showindices(io, axes(a)...)
    print(io, ')')
    toplevel && print(io, " with eltype ", eltype(a))
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

Return an `AbstractArray` that shares structure and has the same type and size as the
argument, but has 1-based indexing. May just return the argument when applicable. Not
exported.

The default implementation uses `OffsetArrays`, but other types should use something more
specific to remove a level of indirection when applicable.

```jldoctest; setup=:(using OffsetArrays)
julia> A = [1 3 5; 2 4 6];

julia> O = OffsetArray(A, 0:1, -1:1)
2×3 OffsetArray(::Matrix{Int64}, 0:1, -1:1) with eltype Int64 with indices 0:1×-1:1:
 1  3  5
 2  4  6

julia> OffsetArrays.no_offset_view(O)[1,1] = -9
-9

julia> A
2×3 Matrix{Int64}:
 -9  3  5
  2  4  6
```
"""
function no_offset_view(A::AbstractArray)
    if Base.has_offset_axes(A)
        OffsetArray(A, map(r->1-first(r), axes(A)))
    else
        A
    end
end

no_offset_view(A::OffsetArray) = no_offset_view(parent(A))


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

end # module
