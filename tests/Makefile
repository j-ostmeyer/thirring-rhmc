#FC = /usr/lib64/mpi/gcc/mvapich2/bin/mpif90
#FC = /usr/lib64/mvapich2/bin/mpif90
MPIFC = mpiifort
FC = ifort
#FC = gfortran
#FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -heap-arrays -g
#FCFLAGS = -g -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1 #-heap-arrays -CB -traceback
#FCFLAGS = -g -O3 -march=native -mtune=native -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1
FCFLAGS = -g -O0
#FCFLAGS = -g -O0 -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1
#FCFLAGS = -O0 -g -DMPI -DNP_X=2 -DNP_Y=2 -DNP_T=2 -CB

MPI=no
NP_X=1
NP_Y=1
NP_T=1

ifeq ($(MPI), yes)
	COMMS_FLAGS = -DMPI -DNP_X=$(NP_X) -DNP_Y=$(NP_Y) -DNP_T=$(NP_T)
	COMMS_LIB = comms.o
	FC = $(MPIFC)
else
	COMMS_LIB = uncomms.o
endif

EXES = test_dslash test_dslashd test_dslash2d test_load test_save test_derivs test_qmrherm test_hamilton test_force test_congrad test_measure

default: $(EXES)

clean:
	rm -f $(EXES) test_halo *.mod *.o

%.o: ../%.f90 Makefile
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) -c -o $@ $<

%.o: ../%.F90 Makefile
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) -c -o $@ $<

test_halo_%: params.o ${COMMS_LIB} test_halo_%.F90
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) -o $@ $^

$(EXES) : % : params.o ${COMMS_LIB} random.o bulk_rhmc_lib.o %.F90
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) $(COMMS_FLAGS) > $@
