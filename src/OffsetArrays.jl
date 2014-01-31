module OffsetArrays

import Base: Array

export OffsetArray

type OffsetArray{T<:Number, N, A<:AbstractArray} <: AbstractArray

    offsets::NTuple{N,Int}
    array::A
    o1::Int
    o2::Int
    o3::Int
    o4::Int
    o5::Int

    function OffsetArray(r::Range1{Int})
        array = Array(T, length(r))
        offs = (1 - minimum(r),)
        new(offs, array, offs[1])
    end

    function OffsetArray(r1::Range1{Int}, r2::Range1{Int})
        dims = (length(r1), length(r2))
        array = Array(T, dims)
        offs = (1 - minimum(r1), 1 - minimum(r2))
        new(offs, array, offs[1], offs[2])
    end

    function OffsetArray(r1::Range1{Int}, r2::Range1{Int}, r3::Range1{Int})
        dims = (length(r1), length(r2), length(r3))
        array = Array(T, dims)
        offs = (1 - minimum(r1), 1 - minimum(r2), 1 - minimum(r3))
        new(offs, array, offs[1], offs[2], offs[3])
    end

    function OffsetArray(r1::Range1{Int}, r2::Range1{Int}, r3::Range1{Int}, r4::Range1{Int})
        dims = (length(r1), length(r2), length(r3), length(r4))
        array = Array(T, dims)
        offs = (1 - minimum(r1), 1 - minimum(r2), 1 - minimum(r3), 1 - minimum(r4))
        new(offs, array, offs[1], offs[2], offs[3], offs[4])
    end

    function OffsetArray(r1::Range1{Int}, r2::Range1{Int}, r3::Range1{Int}, r4::Range1{Int}, r5::Range1{Int})
        dims = (length(r1), length(r2), length(r3), length(r4), length(r5))
        array = Array(T, dims)
        offs = (1 - minimum(r1), 1 - minimum(r2), 1 - minimum(r3), 1 - minimum(r4), 1 - minimum(r5))
        new(offs, array, offs[1], offs[2], offs[3], offs[4], offs[5])
    end

    function OffsetArray(r::Range1{Int}...)
        dims = map((x) -> length(x), r)
        array = Array(T, dims)
        offs = map((x) -> 1 - minimum(x), r)
        new(offs, array)
    end

end

OffsetArray(T, r::Range1{Int}...) = OffsetArray{T, length(r,), Array{T, length(r,)}}(r...)

getindex{T<:Number}(FA::OffsetArray{T,1}, i1::Int) = FA.array[i1+FA.o1]
getindex{T<:Number}(FA::OffsetArray{T,2}, i1::Int, i2::Int) = FA.array[i1+FA.o1, i2+FA.o2]
getindex{T<:Number}(FA::OffsetArray{T,3}, i1::Int, i2::Int, i3::Int) = FA.array[i1+FA.o1, i2+FA.o2, i3+FA.o3]
getindex{T<:Number}(FA::OffsetArray{T,4}, i1::Int, i2::Int, i3::Int, i4::Int) = FA.array[i1+FA.o1, i2+FA.o2, i3+FA.o3, i4+FA.o4]
getindex{T<:Number}(FA::OffsetArray{T,5}, i1::Int, i2::Int, i3::Int, i4::Int, i5::Int) = FA.array[i1+FA.o1, i2+FA.o2, i3+FA.o3, i4+FA.o4, i5+FA.o5]

# a generic not very efficient case    
getindex{T<:Number,N}(FA::OffsetArray{T,N}, I::Int...) = let ind = [I[i] + FA.offsets[i] for i = 1:length(I)]; return FA.array[ind...] end

setindex!{T<:Number}(FA::OffsetArray{T,1}, x, i1::Int) = arrayset(FA.array, convert(T,x), i1+FA.o1)
setindex!{T<:Number}(FA::OffsetArray{T,2}, x, i1::Int, i2::Int) = arrayset(FA.array, convert(T,x), i1+FA.o1, i2+FA.o2)
setindex!{T<:Number}(FA::OffsetArray{T,3}, x, i1::Int, i2::Int, i3::Int) = arrayset(FA.array, convert(T,x), i1+FA.o1, i2+FA.o2, i3+FA.o3)
setindex!{T<:Number}(FA::OffsetArray{T,4}, x, i1::Int, i2::Int, i3::Int, i4::Int) = arrayset(FA.array, convert(T,x), i1+FA.o1, i2+FA.o2, i3+FA.o3, i4+FA.o4)
setindex!{T<:Number}(FA::OffsetArray{T,5}, x, i1::Int, i2::Int, i3::Int, i4::Int, i5::Int) = arrayset(FA.array, convert(T,x), i1+FA.o1, i2+FA.o2, i3+FA.o3, i4+FA.o4, i5+FA.o5)

# a generic not very efficient case    
setindex!{T<:Number,N}(FA::OffsetArray{T,N}, x, I::Int...) = let ind = [I[i] + FA.offsets[i] for i = 1:length(I)]; arrayset(FA.array, convert(T,x), ind...) end

Base.print(a::OffsetArray) = Base.print(a.array)
Base.display(a::OffsetArray) = Base.display(a.array)

Base.size(a::OffsetArray) = arraysize(a.array)
Base.size(a::OffsetArray, d) = arraysize(a.array, d)

end # module
