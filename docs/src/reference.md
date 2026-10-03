# Reference

```@docs
OffsetArray
OffsetVector
OffsetMatrix
OffsetArrays.Origin
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
