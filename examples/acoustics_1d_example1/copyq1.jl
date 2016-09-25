using OffsetArrays     
     
#============================================================================================
function copyq1(meqn::Int,mbc::Int,mx::Int,q1::OffsetArray{Float64},q2::OffsetArray{Float64})
#============================================================================================
     
     # copy the contents of q1 into q2
     
#      implicit double precision (a-h,o-z)
#      dimension q1(meqn,1-mbc:mx+mbc)
#      dimension q2(meqn,1-mbc:mx+mbc)
     
    for i = 1-mbc : mx+mbc
        for m=1:meqn
            q2[m,i] = q1[m,i]
        end
    end

end

