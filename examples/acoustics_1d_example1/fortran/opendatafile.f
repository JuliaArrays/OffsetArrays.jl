
      subroutine opendatafile(iunit, fname)
c     #
c     # Open the file fname and determine how many leading lines are
c     # comments.  Then rewind and skip over the comment lines so that the
c     # file is ready for reading data from.
c     #
c     # All comment lines must start with # in the first column.

      integer iunit, commentlines


c     This allows for a variable name length and obsolves (I think) the
c     length problem that appeared for some compilers (notably ifort)
      character(len=*) :: fname

      character*1 firstchar
      logical foundFile
      
      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
          if (.not. foundFile) then
            print "(2a)",'*** in opendatafile, file not found:', fname 
            stop
          endif
      else
          open(unit=iunit,file=fname,status='old',form='formatted')
          print "(2a)",'Reading data file: ', fname
      endif

      firstchar = '#'
      commentlines = -1
      do commentlines=-1,1000
          if (firstchar .eq. '#') then
              read(iunit,"(a1)") firstchar
          else
              exit
          endif
      enddo
      
      write(6,602) commentlines
  602 format('         first',i2,
     &       ' lines are comments and will be skipped')
     
      return
      end
