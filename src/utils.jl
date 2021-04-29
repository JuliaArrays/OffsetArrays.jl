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

@inline _maybewrapoffset(r::AbstractVector, ax::Tuple{Any}) = _maybewrapoffset(r, ax[1])
@inline _maybewrapoffset(r::AbstractUnitRange{<:Integer}, ::Base.OneTo) = no_offset_view(r)
@inline _maybewrapoffset(r::AbstractVector, ::Base.OneTo) = no_offset_view(r)
@inline function _maybewrapoffset(r::AbstractUnitRange{<:Integer}, ax::AbstractUnitRange)
	of = first(ax) - 1
	IdOffsetRange(_shiftedUnitRange(r, of), of)
end
@inline _maybewrapoffset(r::AbstractVector, ax::AbstractUnitRange) = OffsetArray(r, ax)

_shiftedUnitRange(r::AbstractUnitRange, of) = UnitRange(first(r) - of, last(r) - of)

if VERSION <= v"1.7.0-DEV.1039"
    _contiguousindexingtype(r::AbstractUnitRange{<:Integer}) = UnitRange{Int}(r)
else
    _contiguousindexingtype(r::AbstractUnitRange{<:Integer}) = r
end
