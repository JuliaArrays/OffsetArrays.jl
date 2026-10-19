# Reference

```@docs
OffsetArray
OffsetVector
OffsetMatrix
OffsetArrays.Origin
OffsetArrays.AbstractOffsetArray
OffsetArrays.AbstractOffsetVector
OffsetArrays.AbstractOffsetMatrix
OffsetArrays.offsets
OffsetArrays.parenttype
OffsetArrays.overflow_check
OffsetArrays.unwrap
OffsetArrays.IdOffsetRange
OffsetArrays.no_offset_view
OffsetArrays.AxisConversionStyle
OffsetArrays.center
OffsetArrays.centered
```

# Internals
```@docs
Base.unsafe_wrap(::OffsetArrays.OffsetArrayUnion{T,N}, ::Ptr{T}, ::NTuple{N, OffsetArrays.OffsetAxisKnownLength}) where {T,N}
```
