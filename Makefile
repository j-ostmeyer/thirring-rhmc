FC = ifort
#FC = gfortran
FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -g #-heap-arrays
#FCFLAGS = -O0 -g -heap-arrays

default: bulk_rhmc compile_flags

.PHONY: clean

clean:
	rm -f bulk_rhmc *.o *.mod compile_flags

%.o: %.f90 Makefile
	$(FC) $(FCFLAGS) -c -o $@ $<

%.o: %.F90 Makefile
	$(FC) $(FCFLAGS) -c -o $@ $<

bulk_rhmc: random.o bulk_rhmc_lib.o bulk_rhmc.F90
	$(FC) $(FCFLAGS) -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) > $@
