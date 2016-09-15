# to run this program after it has been made, type
#   flinearad > output
# to see the results after running, type
#   xmgrace output
FC=gfortran
flinearad : linaddmain.o consdiff.o linearad.o riemprob.o upwind.o
	$(FC) -o $@ linaddmain.o consdiff.o linearad.o riemprob.o \
	  upwind.o
clean :
	rm -f *.o flinearad
