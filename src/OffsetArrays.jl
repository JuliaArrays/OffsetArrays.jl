module OffsetArrays

Base.@deprecate_binding (..) Colon()

# FIXME there is no SimIdx anymore
# TODO look into https://github.com/JuliaLang/julia/commit/c318a9f437b9a0c1cf3c12d9c3caac2cd954cbf9#diff-2264bb51acec4e7e2219a3cb1c733651 and fix the code

using Base: SimIdx, Indices, LinearSlow, LinearFast

export OffsetArray

immutable OffsetArray{T,N,AA<:AbstractArray} <: AbstractArray{T,N}
    parent::AA
    offsets::NTuple{N,Int}
end
typealias OffsetVector{T,AA<:AbstractArray} OffsetArray{T,1,AA}

OffsetArray{T,N}(A::AbstractArray{T,N}, offsets::NTuple{N,Int}) = OffsetArray{T,N,typeof(A)}(A, offsets)
OffsetArray{T,N}(A::AbstractArray{T,N}, offsets::Vararg{Int,N}) = OffsetArray(A, offsets)

(::Type{OffsetArray{T,N}}){T,N}(inds::Indices{N}) = OffsetArray{T,N,Array{T,N}}(Array{T,N}(map(Base.dimlength, inds)), map(indsoffset, inds))
(::Type{OffsetArray{T}}){T,N}(inds::Indices{N}) = OffsetArray{T,N}(inds)
OffsetArray{T,N}(::Type{T}, inds::Vararg{UnitRange{Int},N}) = OffsetArray{T,N}(inds)

Base.linearindexing{T<:OffsetArray}(::Type{T}) = Base.linearindexing(parenttype(T))
Base.indicesbehavior{T<:OffsetArray}(::Type{T}) = Base.IndicesUnitRange()
parenttype{T,N,AA}(::Type{OffsetArray{T,N,AA}}) = AA
parenttype(A::OffsetArray) = parenttype(typeof(A))

Base.parent(A::OffsetArray) = A.parent
Base.size(A::OffsetArray) = size(parent(A))
Base.eachindex(::LinearSlow, A::OffsetArray) = CartesianRange(indices(A))
Base.eachindex(::LinearFast, A::OffsetVector) = indices(A, 1)
Base.summary(A::OffsetArray) = string(typeof(A))*" with indices "*string(indices(A))

# Implementations of indices and indices1. Since bounds-checking is
# performance-critical and relies on indices, these are usually worth
# optimizing thoroughly.
@inline Base.indices(A::OffsetArray, d) = 1 <= d <= length(A.offsets) ? (o = A.offsets[d]; (1+o:size(parent(A),d)+o)) : (1:1)
@inline Base.indices(A::OffsetArray) = _indices((), A)  # would rather use ntuple, but see #15276
@inline _indices{T,N}(out::NTuple{N}, A::OffsetArray{T,N}) = out
@inline _indices(out, A::OffsetArray) = (d = length(out)+1; o = A.offsets[d]; _indices((out..., (1+o:size(parent(A),d)+o)), A))
# By optimizing indices1 we can avoid a branch on the dim-check
Base.indices1{T}(A::OffsetArray{T,0}) = 1:1
@inline Base.indices1(A::OffsetArray) = (o = A.offsets[1]; 1+o:size(parent(A),1)+o)

function Base.similar(A::OffsetArray, T::Type, dims::Dims)
    B = similar(parent(A), T, dims)
end
function Base.similar(A::AbstractArray, T::Type, inds::Tuple{Vararg{SimIdx}})
    B = similar(A, T, map(Base.dimlength, inds))
    OffsetArray(B, map(indsoffset, inds))
end

Base.allocate_for(f, A::OffsetArray, shape::SimIdx) = OffsetArray(f(Base.dimlength(shape)), (indsoffset(shape),))
Base.allocate_for(f, A::OffsetArray, shape::Tuple{Vararg{SimIdx}}) = OffsetArray(f(map(Base.dimlength, shape)), map(indsoffset, shape))
Base.promote_indices(a::OffsetArray, b::OffsetArray) = a

Base.reshape(A::AbstractArray, inds::Indices) = OffsetArray(reshape(A, map(Base.dimlength, inds)), map(indsoffset, inds))

@inline function Base.getindex{T,N}(A::OffsetArray{T,N}, I::Vararg{Int,N})
    @boundscheck checkbounds(A, I...)
    @inbounds ret = parent(A)[offset(A.offsets, I)...]
    ret
end
@inline function Base._getindex(::LinearFast, A::OffsetVector, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds ret = parent(A)[offset(A.offsets, (i,))[1]]
    ret
end
@inline function Base._getindex(::LinearFast, A::OffsetArray, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds ret = parent(A)[i]
    ret
end
@inline function Base.setindex!{T,N}(A::OffsetArray{T,N}, val, I::Vararg{Int,N})
    @boundscheck checkbounds(A, I...)
    @inbounds parent(A)[offset(A.offsets, I)...] = val
    val
end
@inline function Base._setindex!(::LinearFast, A::OffsetVector, val, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[offset(A.offsets, (i,))[1]] = val
    val
end
@inline function Base._setindex!(::LinearFast, A::OffsetArray, val, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[i] = val
    val
end

# Computing a shifted index (subtracting the offset)
offset{N}(offsets::NTuple{N,Int}, inds::NTuple{N,Int}) = _offset((), offsets, inds)
_offset(out, ::Tuple{}, ::Tuple{}) = out
@inline _offset(out, offsets, inds) = _offset((out..., inds[1]-offsets[1]), Base.tail(offsets), Base.tail(inds))

indsoffset(r::Range) = first(r) - 1
indsoffset(i::Integer) = 0

end # module
