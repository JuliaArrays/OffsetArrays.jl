module OffsetArrays

#import Base: Array
importall Base

export OffsetArray, ..

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

# a generic but not very efficient case    
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


# as the parser changes colons into ranges.
# for avoiding to write a[Colon()] a whole range .. is introduced
# it is possible to write say a[..] = b[..] 

const (..) = Colon()

getindex{T<:Number}(FA::OffsetArray{T,1}, r1::Union(Range1{Int},Colon)) =
    isa(r1, Colon) ? FA.array[:] : FA.array[r1+FA.o1]

getindex{T<:Number}(FA::OffsetArray{T,2}, r1::Union(Range1{Int},Colon), i2::Int) =
    isa(r1, Colon) ? FA.array[:, i2+FA.o2] : FA.array[r1+FA.o1, i2+FA.o2]

getindex{T<:Number}(FA::OffsetArray{T,3}, r1::Union(Range1{Int},Colon), i2::Int, i3::Int) =
    isa(r1, Colon) ? FA.array[:, i2+FA.o2, i3+FA.o3] : FA.array[r1+FA.o1, i2+FA.o2, i3+FA.o3]

setindex!{T<:Number}(FA::OffsetArray{T,1}, x, r1::Union(Range1{Int},Colon)) = let
    if isa(r1, Colon)
        FA.array[:] = x[:]
    else
        FA.array[r1+FA.o1] = x[r1+FA.o1]
    end
end

setindex!{T<:Number}(FA::OffsetArray{T,2}, x, r1::Union(Range1{Int},Colon), i2::Int) = let
    if isa(r1, Colon)
        FA.array[:,        i2+FA.o2] = x[:,        i2+FA.o2]
    else
        FA.array[r1+FA.o1, i2+FA.o2] = x[r1+FA.o1, i2+FA.o2]
    end
end

setindex!{T<:Number}(FA::OffsetArray{T,3}, x, r1::Union(Range1{Int},Colon), i2::Int, i3::Int) = let
    if isa(r1, Colon)
        FA.array[:,        i2+FA.o2, i3+FA.o3] = x[:,        i2+FA.o2, i3+FA.o3]
    else
        FA.array[r1+FA.o1, i2+FA.o2, i3+FA.o3] = x[r1+FA.o1, i2+FA.o2, i3+FA.o3]
    end
end


end # module
