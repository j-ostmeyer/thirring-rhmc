FC = ifort
#FC = gfortran
#FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -heap-arrays -g
FCFLAGS = -O0 -g -heap-arrays

EXES = test_dslash

default: $(EXES)

clean:
	rm bulk_rhmc_lib.o $(EXES)

bulk_rhmc_lib.o: ../bulk_rhmc_lib.f
	$(FC) $(FCFLAGS) -c -o $@ $^

$(EXES) : % : %.f90 bulk_rhmc_lib.o
	$(FC) $(FCFLAGS) -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) > $@
