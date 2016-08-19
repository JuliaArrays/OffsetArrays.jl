module OffsetArrays

Base.@deprecate_binding (..) Colon()

using Base: Indices, LinearSlow, LinearFast, tail

export OffsetArray

immutable OffsetArray{T,N,AA<:AbstractArray} <: AbstractArray{T,N}
    parent::AA
    offsets::NTuple{N,Int}
end
typealias OffsetVector{T,AA<:AbstractArray} OffsetArray{T,1,AA}

OffsetArray{T,N}(A::AbstractArray{T,N}, offsets::NTuple{N,Int}) =
    OffsetArray{T,N,typeof(A)}(A, offsets)
OffsetArray{T,N}(A::AbstractArray{T,N}, offsets::Vararg{Int,N}) =
    OffsetArray(A, offsets)

(::Type{OffsetArray{T,N}}){T,N}(inds::Indices{N}) =
    OffsetArray{T,N,Array{T,N}}(Array{T,N}(map(length, inds)), map(indexoffset, inds))
(::Type{OffsetArray{T}}){T,N}(inds::Indices{N}) = OffsetArray{T,N}(inds)
OffsetArray{T}(A::AbstractArray{T,0}) = OffsetArray{T,0,typeof(A)}(A, ())
OffsetArray{T,N}(::Type{T}, inds::Vararg{UnitRange{Int},N}) = OffsetArray{T,N}(inds)

# The next two are necessary for ambiguity resolution. Really, the
# second method should not be necessary.
OffsetArray{T}(A::AbstractArray{T,0}, inds::Tuple{}) = OffsetArray{T,0,typeof(A)}(A, ())
OffsetArray{T,N}(A::AbstractArray{T,N}, inds::Tuple{}) = error("this should never be called")
OffsetArray{T,N}(A::AbstractArray{T,N}, inds::NTuple{N,AbstractUnitRange}) =
    OffsetArray(A, map(indexoffset, inds))
OffsetArray{T,N}(A::AbstractArray{T,N}, inds::Vararg{AbstractUnitRange,N}) =
    OffsetArray(A, inds)

Base.linearindexing{T<:OffsetArray}(::Type{T}) = Base.linearindexing(parenttype(T))
parenttype{T,N,AA}(::Type{OffsetArray{T,N,AA}}) = AA
parenttype(A::OffsetArray) = parenttype(typeof(A))

Base.parent(A::OffsetArray) = A.parent

errmsg(A) = error("size not supported for arrays with indices $(indices(A)); see http://docs.julialang.org/en/latest/devdocs/offset-arrays/")
Base.size(A::OffsetArray) = errmsg(A)
Base.size(A::OffsetArray, d) = errmsg(A)
Base.eachindex(::LinearSlow, A::OffsetArray) = CartesianRange(indices(A))
Base.eachindex(::LinearFast, A::OffsetVector) = indices(A, 1)

# Implementations of indices and indices1. Since bounds-checking is
# performance-critical and relies on indices, these are usually worth
# optimizing thoroughly.
@inline Base.indices(A::OffsetArray, d) = 1 <= d <= length(A.offsets) ? indices(parent(A))[d] + A.offsets[d] : (1:1)
@inline Base.indices(A::OffsetArray) = _indices(indices(parent(A)), A.offsets)  # would rather use ntuple, but see #15276
@inline _indices(inds, offsets) = (inds[1]+offsets[1], _indices(tail(inds), tail(offsets))...)
_indices(::Tuple{}, ::Tuple{}) = ()
Base.indices1{T}(A::OffsetArray{T,0}) = 1:1  # we only need to specialize this one

function Base.similar(A::OffsetArray, T::Type, dims::Dims)
    B = similar(parent(A), T, dims)
end
function Base.similar(A::AbstractArray, T::Type, inds::Tuple{UnitRange,Vararg{UnitRange}})
    B = similar(A, T, map(length, inds))
    OffsetArray(B, map(indexoffset, inds))
end

Base.similar(f::Union{Function,DataType}, shape::Tuple{UnitRange,Vararg{UnitRange}}) = OffsetArray(f(map(length, shape)), map(indexoffset, shape))

Base.reshape(A::AbstractArray, inds::Tuple{UnitRange,Vararg{UnitRange}}) = OffsetArray(reshape(A, map(length, inds)), map(indexoffset, inds))

Base.reshape(A::OffsetArray, inds::Tuple{UnitRange,Vararg{UnitRange}}) = OffsetArray(reshape(parent(A), map(length, inds)), map(indexoffset, inds))
# Don't allow bounds-checks to be removed during Julia 0.5
@inline function Base.getindex{T,N}(A::OffsetArray{T,N}, I::Vararg{Int,N})
    checkbounds(A, I...)
    @inbounds ret = parent(A)[offset(A.offsets, I)...]
    ret
end
@inline function Base._getindex(::LinearFast, A::OffsetVector, i::Int)
    checkbounds(A, i)
    @inbounds ret = parent(A)[offset(A.offsets, (i,))[1]]
    ret
end
@inline function Base._getindex(::LinearFast, A::OffsetArray, i::Int)
    checkbounds(A, i)
    @inbounds ret = parent(A)[i]
    ret
end
@inline function Base.setindex!{T,N}(A::OffsetArray{T,N}, val, I::Vararg{Int,N})
    checkbounds(A, I...)
    @inbounds parent(A)[offset(A.offsets, I)...] = val
    val
end
@inline function Base._setindex!(::LinearFast, A::OffsetVector, val, i::Int)
    checkbounds(A, i)
    @inbounds parent(A)[offset(A.offsets, (i,))[1]] = val
    val
end
@inline function Base._setindex!(::LinearFast, A::OffsetArray, val, i::Int)
    checkbounds(A, i)
    @inbounds parent(A)[i] = val
    val
end

### Convenience functions ###

Base.zeros(T::Type, inds::UnitRange...) = fill!(OffsetArray{T}(inds), zero(T))
Base.zeros(inds::UnitRange...) = zeros(Float64, inds...)
Base.ones(T::Type, inds::UnitRange...) = fill!(OffsetArray{T}(inds), one(T))
Base.ones(inds::UnitRange...) = ones(Float64, inds...)
Base.fill(x, inds::Tuple{UnitRange,Vararg{UnitRange}}) = fill!(OffsetArray{typeof(x)}(inds), x)
@inline Base.fill(x, ind1::UnitRange, inds::UnitRange...) = fill(x, (ind1, inds...))

### Low-level utilities ###

# Computing a shifted index (subtracting the offset)
offset{N}(offsets::NTuple{N,Int}, inds::NTuple{N,Int}) = _offset((), offsets, inds)
_offset(out, ::Tuple{}, ::Tuple{}) = out
@inline _offset(out, offsets, inds) = _offset((out..., inds[1]-offsets[1]), Base.tail(offsets), Base.tail(inds))

indexoffset(r::Range) = first(r) - 1
indexoffset(i::Integer) = 0

end # module
