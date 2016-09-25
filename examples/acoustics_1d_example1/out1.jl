using OffsetArrays
#import Base.close

#=========================================================#
function out1(meqn::Int, mbc::Int, mx::Int, xlower::Float64, dx::Float64,
              q::OffsetArray{Float64,2}, t::Float64, iframe::Int,
              aux::OffsetArray{Float64}, maux::Int, outaux::Bool)
#=========================================================#
#     
#     # Output the results for a general system of conservation laws
#     # in 1 dimension
#     
#     # Write the results to the file fort.q<iframe>
#     # Use format required by matlab script  plotclaw1.m
#     # The same format is used by the amrclaw package.
#     # Here it's adapted to output just the single grid.
#     
#     # set outaux = .true. to also output the aux arrays to fort.a<iframe>
#     
#     
##       implicit double precision (a-h,o-z)
##       dimension q(meqn,1-mbc:mx+mbc)
##       dimension aux(maux,1-mbc:mx+mbc)
##       character*10 fname1, fname2, fname3
##       logical outaux

#     
#     # Write the results to the file fort.q<iframe>
#     # Use format required by matlab script  plotclaw1.m
#     # The same format is used by the amrclaw package.  
#     # Here it's adapted to output just the single grid.
#     
#     # first create the file name and open file
#     
##      fname1 = 'fort.qxxxx'
##      fname2 = 'fort.txxxx'
##      fname3 = 'fort.axxxx'
##      nstp = iframe
##      do 55 ipos = 10, 7, -1
##         idigit = mod(nstp,10)
##         fname1(ipos:ipos) = char(ichar('0') + idigit)
##         fname2(ipos:ipos) = char(ichar('0') + idigit)
##         fname3(ipos:ipos) = char(ichar('0') + idigit)
##         nstp = nstp / 10
## 55   continue
##

    fname1 = @sprintf("fort.q%04d", iframe)
    fname2 = @sprintf("fort.t%04d", iframe)
    fname3 = @sprintf("fort.a%04d", iframe)

##      open(unit=50,file=fname1,status='unknown',form='formatted')
##      open(unit=60,file=fname2,status='unknown',form='formatted')

    unit50 = open(fname1, "w")
    unit60 = open(fname2, "w")

    # the following parameters are used in amrclaw where there are
    # multiple grids.  Here they are all set to 1:
    ngrids = 1
    mptr = 1
    level = 1

    # write(50,1001) mptr,level,mx

## 1001 format(i5,'                 grid_number',/,
##     &     i5,'                 AMR_level',/,
##     &     i5,'                 mx')
##

    @printf(unit50,"%5d                 grid_number\n",mptr)
    @printf(unit50,"%5d                 AMR_level\n",level)
    @printf(unit50,"%5d                 mx\n",mx)
    
##      write(50,1002) xlower,dx
## 1002 format(e18.8,'    xlow', /,
##     &     e18.8,'    dx', /)
###     

    @printf(unit50,"%18.8e    xlow\n",xlower)
    @printf(unit50,"%18.8e    dx\n",dx)

##      do 10 i=1,mx
##         do m=1,meqn
###     # exponents with more than 2 digits cause problems reading
###     # into matlab... reset tiny values to zero:
##            if (dabs(q(m,i)) .lt. 1d-99) q(m,i) = 0.d0
##         enddo
###     
##         write(50,1005) (q(m,i), m=1,meqn)
## 1005    format(4e16.8)
###     
## 10   continue
##      write(50,*) ' '
##      write(50,*) ' '
##

    for i=1:mx
        for m=1:meqn
            if abs(q[m,i]) < 1e-99
                q[m,i] = 0.0
            end
            @printf(unit50,"%16.8e",q[m,i])
            # FIXME
            ## if i % 4 == 1 && i > 1
            ##     @printf(unit50,"\n")
            ## end
        end
        @printf(unit50,"\n")
    end
    @printf(unit50," \n \n")


##
##      if (outaux) then 
###     # also output the aux arrays:
##         open(unit=70,file=fname3,status='unknown',form='formatted')
##         write(70,1001) mptr,level,mx
##         write(70,1002) xlower,dx
##         do 110 i=1,mx
##            do m=1,maux
###     # exponents with more than 2 digits cause problems reading
###     # into matlab... reset tiny values to zero:
##               if (dabs(aux(m,i)) .lt. 1d-99) aux(m,i) = 0.d0
##            enddo
###     
##            write(70,1005) (aux(m,i), m=1,maux)
###     
## 110     continue
##         write(70,*) ' '
##         close(unit=70)
##      endif
##
##      write(60,1000) t,meqn,ngrids,maux,1
##
## 1000 format(e26.16,'    time', /,
##     &     i5,'                 meqn'/,
##     &     i5,'                 ngrids'/,
##     &     i5,'                 maux'/,
##     &     i5,'                 ndim'/,/)
###     
##
    Base.close(unit50)
    Base.close(unit60)

end

