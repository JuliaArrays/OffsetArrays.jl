#***********************************************************************
#  Copyright 2006 John A. Trangenstein
#
#  This software is made available for research and instructional use 
#  only. 
#  You may copy and use this software without charge for these 
#  non-commercial purposes, provided that the copyright notice and 
#  associated text is reproduced on all copies.  
#  For all other uses (including distribution of modified versions), 
#  please contact the author at
#    John A. Trangenstein
#    Department of Mathematics
#    Duke University
#    Durham, NC 27708-0320
#    USA
#  or
#    johnt@math.duke.edu
#  
#  This software is made available "as is" without any assurance that it
#  is completely correct, or that it will work for your purposes.  
#  Use the software at your own risk.
#***********************************************************************

module consts

    const roundoff  =  1.0-14
    const small     =  1.0-20
    const huge      =  1.0e300
    const undefind  =  1.0e300
    const lnrndoff  =  14.0e0
    const lnsmall   =  20.0e0

    # problem-specific parameters:
    const jump      =  0.0
    const x_left    = -0.2
    const x_right   =  1.0
    const statelft  =  2.0
    const statergt  =  0.0
    const velocity  =  1.0

    export roundoff, huge, jump, x_left, x_right, statelft, statergt, velocity

end

