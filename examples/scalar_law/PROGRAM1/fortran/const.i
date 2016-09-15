c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
c static char sccsid[] = "%W%\t%G%"
      double precision 
     &  roundoff,small,huge,undefind,lnrndoff,lnsmall
      integer ihuge,inhuge
      common/machine/roundoff,small,huge,undefind,lnrndoff,lnsmall,
     &  ihuge,inhuge
      double precision 
     &  zero,sixtyfourth,thirtysecond,sixteenth,tenth,eighth,sixth,
     &  fourth,third,point4,half,twothird,pt75,onemert2,rt75,one,
     &  fourthir,onept5,two,three,pi,four,seven,nine,ten
      parameter (zero=0.d0)
      parameter (sixtyfourth=0.015625d0)
      parameter (thirtysecond=0.03125d0)
      parameter (sixteenth=0.0625d0)
      parameter (tenth=0.1d0)
      parameter (eighth=0.125d0)
      parameter (sixth=0.16666666666667d0)
      parameter (fourth=.25d0)
      parameter (third=.333333333333333d0)
      parameter (point4=.4d0)
      parameter (half=.5d0)
      parameter (twothird=.66666666666667d0)
      parameter (pt75=.75d0)
      parameter (onemert2=.75688326556578578920d0)
      parameter (rt75=.8660254037844d0)
      parameter (one=1.d0)
      parameter (fourthir=1.33333333333333d0)
      parameter (onept5=1.5d0)
      parameter (two=2.d0)
      parameter (three=3.d0)
      parameter (pi=3.14159265358979323846d0)
      parameter (four=4.d0)
      parameter (seven=7.d0)
      parameter (nine=9.d0)
      parameter (ten=10.d0)
