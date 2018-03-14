FC = /usr/lib64/mpi/gcc/mvapich2/bin/mpif90
#FC = mpiifort
#FC = ifort
#FC = gfortran
#FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -heap-arrays -g
#FCFLAGS = -g -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1 #-heap-arrays -CB -traceback
FCFLAGS = -g -O3 -march=native -mtune=native  -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1
#FCFLAGS = -O0 -g -heap-arrays 

EXES = test_dslash test_dslashd test_dslash2d test_load test_save test_derivs test_qmrherm test_hamilton test_force test_congrad test_measure

default: $(EXES)

clean:
	rm -f random.o bulk_rhmc_lib.o $(EXES) *.mod

%.o: ../%.f90 Makefile
	$(FC) $(FCFLAGS) -c -o $@ $<

%.o: ../%.F90 Makefile
	$(FC) $(FCFLAGS) -c -o $@ $<

$(EXES) : % : %.F90 random.o bulk_rhmc_lib.o
	$(FC) $(FCFLAGS) -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) > $@
