# to run this program after it has been made, type
#   flinearad > output
# to see the results after running, type
#   xmgrace output
FC=gfortran
FCFLAGS=-O4
FLFLAGS=-O4

%.o: %.f
	    $(FC) $(FCFLAGS) -c -o $@ $<

flinearad : linaddmain.o consdiff.o linearad.o riemprob.o upwind.o
	$(FC) $(FLFLAGS) -o $@ linaddmain.o consdiff.o linearad.o riemprob.o \
	  upwind.o

clean :
	rm -f *.o flinearad
