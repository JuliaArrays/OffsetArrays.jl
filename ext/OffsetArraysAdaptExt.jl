module OffsetArraysAdaptExt

using OffsetArrays, Adapt

##
# Adapt allows for automatic conversion of CPU OffsetArrays to GPU OffsetArrays
##
import Adapt
Adapt.adapt_structure(to, O::OffsetArray) = OffsetArrays.parent_call(x -> Adapt.adapt(to, x), O)

end
