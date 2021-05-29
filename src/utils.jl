### Low-level utilities ###

_indexoffset(r::AbstractRange) = first(r) - 1
_indexoffset(i::Integer) = 0
_indexoffset(i::Colon) = 0
_indexlength(r::AbstractRange) = length(r)
_indexlength(i::Integer) = i
_indexlength(i::Colon) = Colon()

_strip_IdOffsetRange(r::IdOffsetRange) = parent(r)
_strip_IdOffsetRange(r) = r

_offset(axparent::AbstractUnitRange, ax::AbstractUnitRange) = first(ax) - first(axparent)
_offset(axparent::AbstractUnitRange, ::Union{Integer, Colon}) = 1 - first(axparent)

_offsets(A::AbstractArray) = map(ax -> first(ax) - 1, axes(A))

_offset_reshape_uncolon(::AbstractArray, I::Tuple{OffsetAxisKnownLength,Vararg{OffsetAxisKnownLength}}) = I
function _offset_reshape_uncolon(A::AbstractArray, I::Tuple{OffsetAxis,Vararg{OffsetAxis}})
    @noinline throw1(I) = throw(DimensionMismatch(string("new dimensions $(I) ",
        "may have at most one omitted dimension specified by `Colon()`")))
    @noinline throw2(A, I) = throw(DimensionMismatch(string("array size $(length(A)) ",
        "must be divisible by the product of the new dimensions $I")))

    pre = _before_colon(I...)
    post = _after_colon(I...)
    _any_colon(post...) && throw1(dims)
    nprev = isempty(pre) ? 1 : mapreduce(_length, *, pre)
    npost = isempty(post) ? 1 : mapreduce(_length, *, post)
    sz, remainder = divrem(length(A), nprev * npost)
    remainder == 0 || throw2(A, dims)

    # preserve the offset information, the default offset beyond `ndims(A)` is 0
    n = length(pre)
    r = Base.OneTo(Int(sz))
    Δ = n < ndims(A) ? _offset(axes(A, n+1), r) : 0
    (pre..., IdOffsetRange(r, -Δ), post...)
end
@inline _any_colon() = false
@inline _any_colon(r::Colon, tail...) = true
@inline _any_colon(r::Any, tail...) = _any_colon(tail...)
@inline _before_colon(r::Any, tail...) = (r, _before_colon(tail...)...)
@inline _before_colon(r::Colon, tail...) = ()
@inline _after_colon(r::Any, tail...) =  _after_colon(tail...)
@inline _after_colon(r::Colon, tail...) = tail
@inline _length(r::AbstractUnitRange) = length(r)
@inline _length(n::Int) = n

"""
    OffsetArrays.AxisConversionStyle(typeof(indices))

`AxisConversionStyle` declares if `indices` should be converted to a single `AbstractUnitRange{Int}`
or to a `Tuple{Vararg{AbstractUnitRange{Int}}}` while flattening custom types into indices.
This method is called after `to_indices(A::Array, axes(A), indices)` to provide
further information in case `to_indices` does not return a `Tuple` of `AbstractUnitRange{Int}`.

Custom index types should extend `AxisConversionStyle` and return either `OffsetArray.SingleRange()`,
which is the default, or `OffsetArray.TupleOfRanges()`. In the former case, the type `T` should
define `Base.convert(::Type{AbstractUnitRange{Int}}, ::T)`, whereas in the latter it should define
`Base.convert(::Type{Tuple{Vararg{AbstractUnitRange{Int}}}}, ::T)`.

An example of the latter is `CartesianIndices`, which is converted to a `Tuple` of
`AbstractUnitRange{Int}` while flattening the indices.

# Example
```jldoctest; setup=:(using OffsetArrays)
julia> struct NTupleOfUnitRanges{N}
           x ::NTuple{N, UnitRange{Int}}
       end

julia> Base.to_indices(A, inds, t::Tuple{NTupleOfUnitRanges{N}}) where {N} = t;

julia> OffsetArrays.AxisConversionStyle(::Type{NTupleOfUnitRanges{N}}) where {N} = OffsetArrays.TupleOfRanges();

julia> Base.convert(::Type{Tuple{Vararg{AbstractUnitRange{Int}}}}, t::NTupleOfUnitRanges) = t.x;

julia> a = zeros(3, 3);

julia> inds = NTupleOfUnitRanges((3:5, 2:4));

julia> oa = OffsetArray(a, inds);

julia> axes(oa, 1) == 3:5
true

julia> axes(oa, 2) == 2:4
true
```
"""
abstract type AxisConversionStyle end
struct SingleRange <: AxisConversionStyle end
struct TupleOfRanges <: AxisConversionStyle end

AxisConversionStyle(::Type) = SingleRange()
AxisConversionStyle(::Type{<:CartesianIndices}) = TupleOfRanges()

_convertTupleAbstractUnitRange(x) = _convertTupleAbstractUnitRange(AxisConversionStyle(typeof(x)), x)
_convertTupleAbstractUnitRange(::SingleRange, x::Int) = (Base.OneTo(x), )
_convertTupleAbstractUnitRange(::SingleRange, x) = (convert(AbstractUnitRange{Int}, x),)
_convertTupleAbstractUnitRange(::TupleOfRanges, x) = convert(Tuple{Vararg{AbstractUnitRange{Int}}}, x)

_toAbstractUnitRanges(t::Tuple) = (_convertTupleAbstractUnitRange(first(t))..., _toAbstractUnitRanges(tail(t))...)
_toAbstractUnitRanges(::Tuple{}) = ()

# ensure that the indices are consistent in the constructor
_checkindices(A::AbstractArray, indices, label) = _checkindices(ndims(A), indices, label)
function _checkindices(N::Integer, indices, label)
    throw_argumenterror(N, indices, label) = throw(ArgumentError(label*" $indices are not compatible with a $(N)D array"))
    N == length(indices) || throw_argumenterror(N, indices, label)
end

@inline _indexedby(r::AbstractVector, ax::Tuple{Any}) = _indexedby(r, ax[1])
@inline _indexedby(r::AbstractUnitRange{<:Integer}, ::Base.OneTo) = no_offset_view(r)
@inline _indexedby(r::AbstractUnitRange{Bool}, ::Base.OneTo) = no_offset_view(r)
@inline _indexedby(r::AbstractVector, ::Base.OneTo) = no_offset_view(r)
@inline function _indexedby(r::AbstractUnitRange{<:Integer}, ax::AbstractUnitRange)
	of = convert(eltype(r), first(ax) - 1)
	IdOffsetRange(_subtractoffset(r, of), of)
end
@inline _indexedby(r::AbstractUnitRange{Bool}, ax::AbstractUnitRange) = OffsetArray(r, ax)
@inline _indexedby(r::AbstractVector, ax::AbstractUnitRange) = OffsetArray(r, ax)

# These functions are equivalent to the broadcasted operation r .- of
# However these ensure that the result is an AbstractRange even if a specific
# broadcasting behavior is not defined for a custom type
_subtractoffset(r::AbstractUnitRange, of) = UnitRange(first(r) - of, last(r) - of)
_subtractoffset(r::AbstractRange, of) = range(first(r) - of, stop = last(r) - of, step = step(r))

if VERSION <= v"1.7.0-DEV.1039"
    _contiguousindexingtype(r::AbstractUnitRange{<:Integer}) = UnitRange{Int}(r)
else
    _contiguousindexingtype(r::AbstractUnitRange{<:Integer}) = r
end

_of_eltype(::Type{T}, M::AbstractArray{T}) where {T} = M
_of_eltype(T, M::AbstractArray) = map(T, M)
