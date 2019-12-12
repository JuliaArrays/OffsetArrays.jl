VERSION < v"0.7.0-beta2.199" && __precompile__()

module OffsetArrays

using Base: Indices, tail, @propagate_inbounds
@static if !isdefined(Base, :IdentityUnitRange)
    const IdentityUnitRange = Base.Slice
else
    using Base: IdentityUnitRange
end

export OffsetArray, OffsetVector

struct OffsetArray{T,N,AA<:AbstractArray} <: AbstractArray{T,N}
    parent::AA
    offsets::NTuple{N,Int}
end
OffsetVector{T,AA<:AbstractArray} = OffsetArray{T,1,AA}

OffsetArray(A::AbstractArray{T,N}, offsets::NTuple{N,Int}) where {T,N} =
    OffsetArray{T,N,typeof(A)}(A, offsets)
OffsetArray(A::AbstractArray{T,N}, offsets::Vararg{Int,N}) where {T,N} =
    OffsetArray(A, offsets)

const ArrayInitializer = Union{UndefInitializer, Missing, Nothing}
OffsetArray{T,N}(init::ArrayInitializer, inds::Indices{N}) where {T,N} =
    OffsetArray{T,N,Array{T,N}}(Array{T,N}(init, map(indexlength, inds)), map(indexoffset, inds))
OffsetArray{T}(init::ArrayInitializer, inds::Indices{N}) where {T,N} = OffsetArray{T,N}(init, inds)
OffsetArray{T,N}(init::ArrayInitializer, inds::Vararg{AbstractUnitRange,N}) where {T,N} = OffsetArray{T,N}(init, inds)
OffsetArray{T}(init::ArrayInitializer, inds::Vararg{AbstractUnitRange,N}) where {T,N} = OffsetArray{T,N}(init, inds)
OffsetArray(A::AbstractArray{T,0}) where {T} = OffsetArray{T,0,typeof(A)}(A, ())

# OffsetVector constructors
OffsetVector(A::AbstractVector, offset) = OffsetArray(A, offset)
OffsetVector{T}(init::ArrayInitializer, inds::AbstractUnitRange) where {T} = OffsetArray{T}(init, inds)

# deprecated constructors
using Base: @deprecate

@deprecate OffsetArray(::Type{T}, inds::Vararg{UnitRange{Int},N}) where {T,N} OffsetArray{T}(undef, inds)
@deprecate OffsetVector(::Type{T}, inds::AbstractUnitRange) where {T} OffsetVector{T}(undef, inds)
@deprecate OffsetArray{T,N}(inds::Indices{N}) where {T,N} OffsetArray{T,N}(undef, inds)
@deprecate OffsetArray{T}(inds::Indices{N}) where {T,N} OffsetArray{T}(undef, inds)
@deprecate OffsetArray{T,N}(inds::Vararg{AbstractUnitRange,N}) where {T,N} OffsetArray{T,N}(undef, inds)
@deprecate OffsetArray{T}(inds::Vararg{AbstractUnitRange,N}) where {T,N} OffsetArray{T}(undef, inds)
@deprecate OffsetVector{T}(inds::AbstractUnitRange) where {T} OffsetVector{T}(undef, inds)

# The next two are necessary for ambiguity resolution. Really, the
# second method should not be necessary.
OffsetArray(A::AbstractArray{T,0}, inds::Tuple{}) where {T} = OffsetArray{T,0,typeof(A)}(A, ())
OffsetArray(A::AbstractArray{T,N}, inds::Tuple{}) where {T,N} = error("this should never be called")
function OffsetArray(A::AbstractArray{T,N}, inds::NTuple{N,AbstractUnitRange}) where {T,N}
    lA = map(indexlength, axes(A))
    lI = map(indexlength, inds)
    lA == lI || throw(DimensionMismatch("supplied axes do not agree with the size of the array (got size $lA for the array and $lI for the indices"))
    OffsetArray(A, map(indexoffset, inds))
end
OffsetArray(A::AbstractArray{T,N}, inds::Vararg{AbstractUnitRange,N}) where {T,N} =
    OffsetArray(A, inds)

# avoid a level of indirection when nesting OffsetArrays
function OffsetArray(A::OffsetArray, inds::NTuple{N,AbstractUnitRange}) where {N}
    OffsetArray(parent(A), inds)
end
OffsetArray(A::OffsetArray{T,0}, inds::Tuple{}) where {T} = OffsetArray{T,0,typeof(A)}(parent(A), ())
OffsetArray(A::OffsetArray{T,N}, inds::Tuple{}) where {T,N} = error("this should never be called")

Base.IndexStyle(::Type{OA}) where {OA<:OffsetArray} = IndexStyle(parenttype(OA))
parenttype(::Type{OffsetArray{T,N,AA}}) where {T,N,AA} = AA
parenttype(A::OffsetArray) = parenttype(typeof(A))

Base.parent(A::OffsetArray) = A.parent

Base.eachindex(::IndexCartesian, A::OffsetArray) = CartesianIndices(axes(A))
Base.eachindex(::IndexLinear, A::OffsetVector)   = axes(A, 1)

Base.size(A::OffsetArray) = size(parent(A))
Base.size(A::OffsetArray, d) = size(parent(A), d)

# Implementations of axes and indices1. Since bounds-checking is
# performance-critical and relies on axes, these are usually worth
# optimizing thoroughly.
@inline Base.axes(A::OffsetArray, d) =
    1 <= d <= length(A.offsets) ? _slice(axes(parent(A))[d], A.offsets[d]) : (1:1)
@inline Base.axes(A::OffsetArray) =
    _axes(axes(parent(A)), A.offsets)  # would rather use ntuple, but see #15276
@inline _axes(inds, offsets) =
    (_slice(inds[1], offsets[1]), _axes(tail(inds), tail(offsets))...)
_axes(::Tuple{}, ::Tuple{}) = ()
Base.axes1(A::OffsetArray{T,0}) where {T} = 1:1  # we only need to specialize this one

# Avoid the kw-arg on the range(r+x, length=length(r)) call in r .+ x
@inline _slice(r, x) = IdentityUnitRange(Base._range(first(r) + x, nothing, nothing, length(r)))

const OffsetAxis = Union{Integer, UnitRange, Base.OneTo, IdentityUnitRange, Colon}
function Base.similar(A::OffsetArray, ::Type{T}, dims::Dims) where T
    B = similar(parent(A), T, dims)
end
function Base.similar(A::AbstractArray, ::Type{T}, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}}) where T
    B = similar(A, T, map(indexlength, inds))
    OffsetArray(B, map(indexoffset, inds))
end

Base.reshape(A::AbstractArray, inds::OffsetAxis...) = reshape(A, inds)
Base.reshape(A::AbstractArray, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}}) =
    OffsetArray(reshape(A, map(indexlength, inds)), map(indexoffset, inds))

# Reshaping OffsetArrays can "pop" the original OffsetArray wrapper and return
# an OffsetArray(reshape(...)) instead of an OffsetArray(reshape(OffsetArray(...)))
Base.reshape(A::OffsetArray, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}}) =
    OffsetArray(reshape(parent(A), map(indexlength, inds)), map(indexoffset, inds))
# And for non-offset axes, we can just return a reshape of the parent directly
Base.reshape(A::OffsetArray, inds::Tuple{Union{Integer,Base.OneTo},Vararg{Union{Integer,Base.OneTo}}}) = reshape(parent(A), inds)
Base.reshape(A::OffsetArray, inds::Dims) = reshape(parent(A), inds)
Base.reshape(A::OffsetArray, ::Colon) = A
Base.reshape(A::OffsetArray, inds::Union{Int,Colon}...) = reshape(parent(A), inds)
Base.reshape(A::OffsetArray, inds::Tuple{Vararg{Union{Int,Colon}}}) = reshape(parent(A), inds)

Base.similar(::Type{T}, shape::Tuple{OffsetAxis,Vararg{OffsetAxis}}) where {T<:AbstractArray} =
    OffsetArray(T(undef, map(indexlength, shape)), map(indexoffset, shape))

Base.fill(v, inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {N} =
    fill!(OffsetArray(Array{typeof(v), N}(undef, map(indexlength, inds)), map(indexoffset, inds)), v)
Base.zeros(::Type{T}, inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {T, N} =
    fill!(OffsetArray(Array{T, N}(undef, map(indexlength, inds)), map(indexoffset, inds)), zero(T))
Base.ones(::Type{T}, inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {T, N} =
    fill!(OffsetArray(Array{T, N}(undef, map(indexlength, inds)), map(indexoffset, inds)), one(T))
Base.trues(inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {N} =
    fill!(OffsetArray(BitArray{N}(undef, map(indexlength, inds)), map(indexoffset, inds)), true)
Base.falses(inds::NTuple{N, Union{Integer, AbstractUnitRange}}) where {N} =
    fill!(OffsetArray(BitArray{N}(undef, map(indexlength, inds)), map(indexoffset, inds)), false)

@inline @propagate_inbounds function Base.getindex(A::OffsetArray{T,N}, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds ret = parent(A)[offset(A.offsets, I)...]
    ret
end
@inline @propagate_inbounds function Base.getindex(A::OffsetVector, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds ret = parent(A)[offset(A.offsets, (i,))[1]]
    ret
end
@inline @propagate_inbounds function Base.getindex(A::OffsetArray, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds ret = parent(A)[i]
    ret
end
@inline @propagate_inbounds function Base.setindex!(A::OffsetArray{T,N}, val, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds parent(A)[offset(A.offsets, I)...] = val
    val
end
@inline @propagate_inbounds function Base.setindex!(A::OffsetVector, val, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[offset(A.offsets, (i,))[1]] = val
    val
end
@inline @propagate_inbounds function Base.setindex!(A::OffsetArray, val, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[i] = val
    val
end

# For fast broadcasting: ref https://discourse.julialang.org/t/why-is-there-a-performance-hit-on-broadcasting-with-offsetarrays/32194
Base.dataids(A::OffsetArray) = Base.dataids(parent(A))

### Special handling for AbstractRange

const OffsetRange{T} = OffsetArray{T,1,<:AbstractRange{T}}
const IIUR = IdentityUnitRange{S} where S<:AbstractUnitRange{T} where T<:Integer

Base.step(a::OffsetRange) = step(parent(a))

Base.getindex(a::OffsetRange, r::OffsetRange) = OffsetArray(a[parent(r)], r.offsets)
Base.getindex(a::OffsetRange, r::AbstractRange) = a.parent[r .- a.offsets[1]]
Base.getindex(a::AbstractRange, r::OffsetRange) = OffsetArray(a[parent(r)], r.offsets)

@inline @propagate_inbounds Base.getindex(r::UnitRange, s::IIUR) =
    OffsetArray(r[s.indices], s)

@inline @propagate_inbounds Base.getindex(r::StepRange, s::IIUR) =
    OffsetArray(r[s.indices], s)

@inline @propagate_inbounds Base.getindex(r::StepRangeLen{T,<:Base.TwicePrecision,<:Base.TwicePrecision}, s::IIUR) where T =
    OffsetArray(r[s.indices], s)
@inline @propagate_inbounds Base.getindex(r::StepRangeLen{T}, s::IIUR) where {T} =
    OffsetArray(r[s.indices], s)

@inline @propagate_inbounds Base.getindex(r::LinRange, s::IIUR) =
    OffsetArray(r[s.indices], s)

function Base.show(io::IO, r::OffsetRange)
    show(io, r.parent)
    o = r.offsets[1]
    print(io, " with indices ", o+1:o+length(r))
end
Base.show(io::IO, ::MIME"text/plain", r::OffsetRange) = show(io, r)

### Convenience functions ###

Base.fill(x, inds::Tuple{UnitRange,Vararg{UnitRange}}) =
    fill!(OffsetArray{typeof(x)}(undef, inds), x)
@inline Base.fill(x, ind1::UnitRange, inds::UnitRange...) = fill(x, (ind1, inds...))


### Some mutating functions defined only for OffsetVector ###

Base.resize!(A::OffsetVector, nl::Integer) = (resize!(A.parent, nl); A)
Base.push!(A::OffsetVector, x...) = (push!(A.parent, x...); A)
Base.pop!(A::OffsetVector) = pop!(A.parent)
Base.empty!(A::OffsetVector) = (empty!(A.parent); A)

### Low-level utilities ###

# Computing a shifted index (subtracting the offset)
@inline offset(offsets::NTuple{N,Int}, inds::NTuple{N,Int}) where {N} =
    (inds[1]-offsets[1], offset(Base.tail(offsets), Base.tail(inds))...)
offset(::Tuple{}, ::Tuple{}) = ()

# Support trailing 1s
@inline offset(offsets::Tuple{Vararg{Int}}, inds::Tuple{Vararg{Int}}) =
    (offset(offsets, Base.front(inds))..., inds[end])
offset(offsets::Tuple{Vararg{Int}}, inds::Tuple{}) = error("inds cannot be shorter than offsets")

indexoffset(r::AbstractRange) = first(r) - 1
indexoffset(i::Integer) = 0
indexoffset(i::Colon) = 0
indexlength(r::AbstractRange) = length(r)
indexlength(i::Integer) = i
indexlength(i::Colon) = Colon()

@eval @deprecate $(Symbol("@unsafe")) $(Symbol("@inbounds"))

function Base.showarg(io::IO, a::OffsetArray, toplevel)
    print(io, "OffsetArray(")
    Base.showarg(io, parent(a), false)
    if ndims(a) > 0
        print(io, ", ")
        printindices(io, axes(a)...)
    end
    print(io, ')')
    toplevel && print(io, " with eltype ", eltype(a))
end
printindices(io::IO, ind1, inds...) =
    (print(io, _unslice(ind1), ", "); printindices(io, inds...))
printindices(io::IO, ind1) = print(io, _unslice(ind1))
_unslice(x) = x
_unslice(x::IdentityUnitRange) = x.indices

"""
    no_offset_view(A)

Return an `AbstractArray` that shares structure and has the same type and size as the
argument, but has 1-based indexing. May just return the argument when applicable. Not
exported.

The default implementation uses `OffsetArrays`, but other types should use something more
specific to remove a level of indirection when applicable.

```jldoctest
julia> O = OffsetArray(A, 0:1, -1:1)
OffsetArray(::Array{Int64,2}, 0:1, -1:1) with eltype Int64 with indices 0:1×-1:1:
 1  3  5
 2  4  6

julia> OffsetArrays.no_offset_view(O)[1,1] = -9
-9

julia> A
2×3 Array{Int64,2}:
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


end # module
