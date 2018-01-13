__precompile__()

module OffsetArrays

using Base: Indices, tail
using Compat
using Compat: axes, CartesianIndices

export OffsetArray, OffsetVector, @unsafe

# TODO: just use .+
# See https://github.com/JuliaLang/julia/pull/22932#issuecomment-330711997
if VERSION < v"0.7.0-DEV.1759"
    plus(r::AbstractUnitRange, i::Integer) = r + i
else
    plus(r::AbstractUnitRange, i::Integer) = broadcast(+, r, i)
end

struct OffsetArray{T,N,AA<:AbstractArray} <: AbstractArray{T,N}
    parent::AA
    offsets::NTuple{N,Int}
end
OffsetVector{T,AA<:AbstractArray} = OffsetArray{T,1,AA}

OffsetArray(A::AbstractArray{T,N}, offsets::NTuple{N,Int}) where {T,N} =
    OffsetArray{T,N,typeof(A)}(A, offsets)
OffsetArray(A::AbstractArray{T,N}, offsets::Vararg{Int,N}) where {T,N} =
    OffsetArray(A, offsets)

OffsetArray{T,N}(inds::Indices{N}) where {T,N} =
    OffsetArray{T,N,Array{T,N}}(Array{T,N}(uninitialized, map(length, inds)), map(indexoffset, inds))
OffsetArray{T}(inds::Indices{N}) where {T,N} = OffsetArray{T,N}(inds)
OffsetArray{T,N}(inds::Vararg{AbstractUnitRange,N}) where {T,N} = OffsetArray{T,N}(inds)
OffsetArray{T}(inds::Vararg{AbstractUnitRange,N}) where {T,N} = OffsetArray{T,N}(inds)
OffsetArray(A::AbstractArray{T,0}) where {T} = OffsetArray{T,0,typeof(A)}(A, ())

# OffsetVector constructors
OffsetVector(A::AbstractVector, offset) = OffsetArray(A, offset)
OffsetVector{T}(inds::AbstractUnitRange) where {T} = OffsetArray{T}(inds)

# https://github.com/JuliaLang/julia/pull/19989
Base.@deprecate OffsetArray(::Type{T}, inds::Vararg{UnitRange{Int},N}) where {T,N} OffsetArray{T}(inds)
Base.@deprecate OffsetVector(::Type{T}, inds::AbstractUnitRange) where {T} OffsetVector{T}(inds)

# The next two are necessary for ambiguity resolution. Really, the
# second method should not be necessary.
OffsetArray(A::AbstractArray{T,0}, inds::Tuple{}) where {T} = OffsetArray{T,0,typeof(A)}(A, ())
OffsetArray(A::AbstractArray{T,N}, inds::Tuple{}) where {T,N} = error("this should never be called")
function OffsetArray(A::AbstractArray{T,N}, inds::NTuple{N,AbstractUnitRange}) where {T,N}
    lA = map(length, axes(A))
    lI = map(length, inds)
    lA == lI || throw(DimensionMismatch("supplied axes do not agree with the size of the array (got size $lA for the array and $lI for the indices"))
    OffsetArray(A, map(indexoffset, inds))
end
OffsetArray(A::AbstractArray{T,N}, inds::Vararg{AbstractUnitRange,N}) where {T,N} =
    OffsetArray(A, inds)

Base.IndexStyle(::Type{OA}) where {OA<:OffsetArray} = IndexStyle(parenttype(OA))
parenttype(::Type{OffsetArray{T,N,AA}}) where {T,N,AA} = AA
parenttype(A::OffsetArray) = parenttype(typeof(A))

Base.parent(A::OffsetArray) = A.parent

errmsg(A) = error("size not supported for arrays with axes $(axes(A)); see http://docs.julialang.org/en/latest/devdocs/offset-arrays/")
Base.size(A::OffsetArray) = errmsg(A)
Base.size(A::OffsetArray, d) = errmsg(A)
Base.eachindex(::IndexCartesian, A::OffsetArray) = CartesianIndices(axes(A))
Base.eachindex(::IndexLinear, A::OffsetVector)   = axes(A, 1)

# Implementations of axes and indices1. Since bounds-checking is
# performance-critical and relies on axes, these are usually worth
# optimizing thoroughly.
@inline Compat.axes(A::OffsetArray, d) =
    1 <= d <= length(A.offsets) ? plus(axes(parent(A))[d], A.offsets[d]) : (1:1)
@inline Compat.axes(A::OffsetArray) =
    _axes(axes(parent(A)), A.offsets)  # would rather use ntuple, but see #15276
@inline _axes(inds, offsets) =
    (plus(inds[1], offsets[1]), _axes(tail(inds), tail(offsets))...)
_axes(::Tuple{}, ::Tuple{}) = ()
Base.indices1(A::OffsetArray{T,0}) where {T} = 1:1  # we only need to specialize this one

function Base.similar(A::OffsetArray, ::Type{T}, dims::Dims) where T
    B = similar(parent(A), T, dims)
end
function Base.similar(A::AbstractArray, ::Type{T}, inds::Tuple{UnitRange,Vararg{UnitRange}}) where T
    B = similar(A, T, map(length, inds))
    OffsetArray(B, map(indexoffset, inds))
end

Base.similar(f::Function, shape::Tuple{UnitRange,Vararg{UnitRange}}) =
    OffsetArray(f(map(length, shape)), map(indexoffset, shape))
Base.similar(::Type{T}, shape::Tuple{UnitRange,Vararg{UnitRange}}) where {T<:OffsetArray} =
    OffsetArray(T(map(length, shape)), map(indexoffset, shape))
Base.similar(::Type{T}, shape::Tuple{UnitRange,Vararg{UnitRange}}) where {T<:Array} =
    OffsetArray(T(uninitialized, map(length, shape)), map(indexoffset, shape))
Base.similar(::Type{T}, shape::Tuple{UnitRange,Vararg{UnitRange}}) where {T<:BitArray} =
    OffsetArray(T(uninitialized, map(length, shape)), map(indexoffset, shape))

Base.reshape(A::AbstractArray, inds::Tuple{UnitRange,Vararg{UnitRange}}) =
    OffsetArray(reshape(A, map(length, inds)), map(indexoffset, inds))

Base.reshape(A::OffsetArray, inds::Tuple{UnitRange,Vararg{UnitRange}}) =
    OffsetArray(reshape(parent(A), map(length, inds)), map(indexoffset, inds))

function Base.reshape(A::OffsetArray, inds::Tuple{UnitRange,Vararg{Union{UnitRange,Int,Base.OneTo}}})
    throw(ArgumentError("reshape must supply UnitRange axes, got $(typeof(inds)).\n       Note that reshape(A, Val{N}) is not supported for OffsetArrays."))
end

# Don't allow bounds-checks to be removed during Julia 0.5
@inline function Base.getindex(A::OffsetArray{T,N}, I::Vararg{Int,N}) where {T,N}
    checkbounds(A, I...)
    @inbounds ret = parent(A)[offset(A.offsets, I)...]
    ret
end
@inline function Base.getindex(A::OffsetVector, i::Int)
    checkbounds(A, i)
    @inbounds ret = parent(A)[offset(A.offsets, (i,))[1]]
    ret
end
@inline function Base.getindex(A::OffsetArray, i::Int)
    checkbounds(A, i)
    @inbounds ret = parent(A)[i]
    ret
end
@inline function Base.setindex!(A::OffsetArray{T,N}, val, I::Vararg{Int,N}) where {T,N}
    checkbounds(A, I...)
    @inbounds parent(A)[offset(A.offsets, I)...] = val
    val
end
@inline function Base.setindex!(A::OffsetVector, val, i::Int)
    checkbounds(A, i)
    @inbounds parent(A)[offset(A.offsets, (i,))[1]] = val
    val
end
@inline function Base.setindex!(A::OffsetArray, val, i::Int)
    checkbounds(A, i)
    @inbounds parent(A)[i] = val
    val
end

### Convenience functions ###

Base.fill(x, inds::Tuple{UnitRange,Vararg{UnitRange}}) =
    fill!(OffsetArray{typeof(x)}(inds), x)
@inline Base.fill(x, ind1::UnitRange, inds::UnitRange...) = fill(x, (ind1, inds...))

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

macro unsafe(ex)
    esc(unsafe(ex))
end
unsafe(ex) = ex
function unsafe(ex::Expr)
    if ex.head ∈ (:+=, :-=, :*=, :/=)
        ex = Expr(:(=), ex.args[1], Expr(:call, Symbol(string(ex.head)[1]), ex.args...))
    end
    if ex.head == :(=)
        a = ex.args[1]
        if isa(a, Expr) && (a::Expr).head == :ref
            # setindex!
            newargs = map(unsafe, ex.args[2:end])
            @assert length(newargs) == 1
            return Expr(:call, :(OffsetArrays.unsafe_setindex!), (a::Expr).args[1], newargs[1], (a::Expr).args[2:end]...)
        end
    end
    newargs = map(unsafe, ex.args)
    if ex.head == :ref
        # getindex
        return Expr(:call, :(OffsetArrays.unsafe_getindex), newargs...)
    end
    Expr(ex.head, newargs...)
end

@inline unsafe_getindex(a::AbstractArray, I...) = (@inbounds ret = a[I...]; ret)
@inline unsafe_setindex!(a::AbstractArray, val, I...) = (@inbounds a[I...] = val; val)

# Linear indexing
@inline unsafe_getindex(a::OffsetArray, i::Int) = _unsafe_getindex(IndexStyle(a), a, i)
@inline unsafe_setindex!(a::OffsetArray, val, i::Int) = _unsafe_setindex!(IndexStyle(a), a, val, i)
for T in (IndexLinear, IndexCartesian)  # ambiguity-resolution requires specificity for both
    @eval begin
        @inline function _unsafe_getindex(::$T, a::OffsetVector, i::Int)
            @inbounds ret = parent(a)[offset(a.offsets, (i,))[1]]
            ret
        end
        @inline function _unsafe_setindex!(::$T, a::OffsetVector, val, i::Int)
            @inbounds parent(a)[offset(a.offsets, (i,))[1]] = val
            val
        end
    end
end
@inline function _unsafe_getindex(::IndexLinear, a::OffsetArray, i::Int)
    @inbounds ret = parent(a)[i]
    ret
end
@inline _unsafe_getindex(::IndexCartesian, a::OffsetArray, i::Int) =
    unsafe_getindex(a, CartesianIndices(axes(a))[i])
@inline function _unsafe_setindex!(::IndexLinear, a::OffsetArray, val, i::Int)
    @inbounds parent(a)[i] = val
    val
end
@inline _unsafe_setindex!(::IndexCartesian, a::OffsetArray, val, i::Int) =
    unsafe_setindex!(a, val, CartesianIndices(axes(a))[i]...)

@inline unsafe_getindex(a::OffsetArray, I::Int...) = unsafe_getindex(parent(a), offset(a.offsets, I)...)
@inline unsafe_setindex!(a::OffsetArray, val, I::Int...) = unsafe_setindex!(parent(a), val, offset(a.offsets, I)...)
@inline unsafe_getindex(a::OffsetArray, I...) = unsafe_getindex(a, Base.IteratorsMD.flatten(I)...)
@inline unsafe_setindex!(a::OffsetArray, val, I...) = unsafe_setindex!(a, val, Base.IteratorsMD.flatten(I)...)

# Indexing a SubArray which has OffsetArray axes
OffsetSubArray{T,N,P,I<:Tuple{OffsetArray,Vararg{OffsetArray}}} = SubArray{T,N,P,I,false}
@inline function unsafe_getindex(a::OffsetSubArray{T,N}, I::Vararg{Int,N}) where {T,N}
    J = map(unsafe_getindex, a.indexes, I)
    unsafe_getindex(parent(a), J...)
end
@inline function unsafe_setindex!(a::OffsetSubArray{T,N}, val, I::Vararg{Int,N}) where {T,N}
    J = map(unsafe_getindex, a.indexes, I)
    unsafe_setindex!(parent(a), val, J...)
end
@inline unsafe_getindex(a::OffsetSubArray, I::Union{Integer,CartesianIndex}...) = unsafe_getindex(a, Base.IteratorsMD.flatten(I)...)
@inline unsafe_setindex!(a::OffsetSubArray, val, I::Union{Integer,CartesianIndex}...) = unsafe_setindex!(a, val, Base.IteratorsMD.flatten(I)...)

if VERSION >= v"0.7.0-DEV.1790"
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
        (print(io, ind1, ", "); printindices(io, inds...))
    printindices(io::IO, ind1) = print(io, ind1)
end

end # module
