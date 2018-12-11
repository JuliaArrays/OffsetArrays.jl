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

const OffsetAxis = Union{Integer, UnitRange, Base.OneTo, IdentityUnitRange}
function Base.similar(A::OffsetArray, ::Type{T}, dims::Dims) where T
    B = similar(parent(A), T, dims)
end
function Base.similar(A::AbstractArray, ::Type{T}, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}}) where T
    B = similar(A, T, map(indexlength, inds))
    OffsetArray(B, map(indexoffset, inds))
end

Base.reshape(A::AbstractArray, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}}) =
    OffsetArray(reshape(A, map(indexlength, inds)), map(indexoffset, inds))

# Reshaping OffsetArrays can "pop" the original OffsetArray wrapper and return
# an OffsetArray(reshape(...)) instead of an OffsetArray(reshape(OffsetArray(...)))
Base.reshape(A::OffsetArray, inds::Tuple{OffsetAxis,Vararg{OffsetAxis}}) =
    OffsetArray(reshape(parent(A), map(indexlength, inds)), map(indexoffset, inds))
# And for non-offset axes, we can just return a reshape of the parent directly
Base.reshape(A::OffsetArray, inds::Tuple{Union{Integer,Base.OneTo},Vararg{Union{Integer,Base.OneTo}}}) = reshape(parent(A), inds)
Base.reshape(A::OffsetArray, inds::Dims) = reshape(parent(A), inds)

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

### Convenience functions ###

Base.fill(x, inds::Tuple{UnitRange,Vararg{UnitRange}}) =
    fill!(OffsetArray{typeof(x)}(undef, inds), x)
@inline Base.fill(x, ind1::UnitRange, inds::UnitRange...) = fill(x, (ind1, inds...))

Base.resize!(A::OffsetVector, nl::Integer) = (resize!(A.parent, nl); A)

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
indexlength(r::AbstractRange) = length(r)
indexlength(i::Integer) = i

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

end # module
