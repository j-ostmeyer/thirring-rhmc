FC = ifort
#FC = gfortran
#FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -heap-arrays -g
FCFLAGS = -g -ipo -O3 -no-prec-div -fp-model fast=2 -xHost #-heap-arrays -CB -traceback
#FCFLAGS = -g -O3 -march=native -mtune=native 
#FCFLAGS = -O0 -g -heap-arrays 

EXES = test_dslash test_dslashd test_dslash2d test_load test_save test_derivs test_qmrherm test_hamilton test_force test_congrad test_measure

default: $(EXES)

clean:
	rm bulk_rhmc_lib.o random.o $(EXES) *.mod

%.o: ../%.f90
	$(FC) $(FCFLAGS) -c -o $@ $^

$(EXES) : % : %.f90 bulk_rhmc_lib.o random.o
	$(FC) $(FCFLAGS) -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) > $@
