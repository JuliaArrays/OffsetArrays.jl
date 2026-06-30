module OffsetArraysAdaptExt

using OffsetArraysCore, Adapt

##
# Adapt allows for automatic conversion of CPU OffsetArrays to GPU OffsetArrays
##
import Adapt
Adapt.adapt_structure(to, O::OffsetArray) = OffsetArraysCore.parent_call(x -> Adapt.adapt(to, x), O)

@static if isdefined(Adapt, :parent_type)
    # To support Adapt 3.0 which doesn't have parent_type defined
    Adapt.parent_type(::Type{OffsetArray{T,N,AA}}) where {T,N,AA} = AA
    Adapt.unwrap_type(W::Type{<:OffsetArray}) = unwrap_type(parent_type(W))
end

end
