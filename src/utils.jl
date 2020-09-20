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

_expandCartesianIndices(inds::Tuple{<:CartesianIndices, Vararg{Any}}) = (convert(Tuple{Vararg{AbstractUnitRange{Int}}}, inds[1])..., _expandCartesianIndices(Base.tail(inds))...)
_expandCartesianIndices(inds::Tuple{Any,Vararg{Any}}) = (inds[1], _expandCartesianIndices(Base.tail(inds))...)
_expandCartesianIndices(::Tuple{}) = ()
