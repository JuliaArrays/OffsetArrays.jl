### Low-level utilities ###

_indexoffset(r::AbstractRange) = first(r) - 1
_indexoffset(i::Integer) = 0
_indexoffset(i::Colon) = 0
_indexlength(r::AbstractRange) = length(r)
_indexlength(i::Integer) = i
_indexlength(i::Colon) = Colon()

_offset(axparent::AbstractUnitRange, ax::AbstractUnitRange) = first(ax) - first(axparent)
_offset(axparent::AbstractUnitRange, ax::CartesianIndices) = _offset(axparent, first(ax.indices))
_offset(axparent::AbstractUnitRange, ax::Integer) = 1 - first(axparent)

_uncolonindices(A::AbstractArray{<:Any,N}, inds::NTuple{N,Any}) where {N} = _uncolonindices(axes(A), inds)
_uncolonindices(ax::Tuple, inds::Tuple) = (first(inds), _uncolonindices(tail(ax), tail(inds))...)
_uncolonindices(ax::Tuple, inds::Tuple{Colon, Vararg{Any}}) = (first(ax), _uncolonindices(tail(ax), tail(inds))...)
_uncolonindices(::Tuple{}, ::Tuple{}) = ()

# Specify offsets using CartesianIndices (issue #71)
# Support a mix of AbstractUnitRanges and CartesianIndices{1}
# Extract the range r from CartesianIndices((r,))
function _stripCartesianIndices(inds::Tuple{CartesianIndices{1},Vararg{Any}})
    I = first(inds)
    Ir = convert(Tuple{AbstractUnitRange{Int}}, I) |> first
    (Ir, _stripCartesianIndices(tail(inds))...)
end
_stripCartesianIndices(inds::Tuple)= (first(inds), _stripCartesianIndices(tail(inds))...)
_stripCartesianIndices(::Tuple{}) = ()

# Support an arbitrary CartesianIndices alongside colons and ranges
# The total number of indices should equal ndims(arr)
# We split the CartesianIndices{N} into N CartesianIndices{1} indices to facilitate dispatch
_splitCartesianIndices(c::CartesianIndices{0}) = ()
function _splitCartesianIndices(c::CartesianIndices)
    c1, ct = Base.IteratorsMD.split(c, Val(1))
    (c1, _splitCartesianIndices(ct)...)
end
function _splitCartesianIndices(t::Tuple{CartesianIndices, Vararg{Any}})
    (_splitCartesianIndices(first(t))..., _splitCartesianIndices(tail(t))...)
end
function _splitCartesianIndices(t::Tuple)
    (first(t), _splitCartesianIndices(tail(t))...)
end
_splitCartesianIndices(::Tuple{}) = ()
