FC = ifort
#FC = gfortran
FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -heap-arrays -g
#FCFLAGS = -O0 -g -heap-arrays

default: bulk_rhmc compile_flags

clean:
	rm bulk_rhmc *.o compile_flags

bulk_rhmc_lib.o: bulk_rhmc_lib.f
	$(FC) $(FCFLAGS) -c -o $@ $^

bulk_rhmc: bulk_rhmc_lib.o bulk_rhmc.f90
	$(FC) $(FCFLAGS) -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) > $@
