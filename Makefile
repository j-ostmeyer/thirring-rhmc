#FC = /usr/lib64/mpi/gcc/mvapich2/bin/mpif90
#FC = /usr/lib64/mvapich2/bin/mpif90
MPIFC = mpiifort
FC = ifort
#FC = gfortran
#FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -heap-arrays -g
#FCFLAGS = -g -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1 #-heap-arrays -CB -traceback
#FCFLAGS = -g -O3 -march=native -mtune=native -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1
FCFLAGS = -g -O0 -CB -heap-arrays -warn all

MPI=yes
NP_X=1
NP_Y=1
NP_T=1

SITE_RANDOM=yes

ifeq ($(MPI), yes)
	COMMS_FLAGS = -DMPI -DNP_X=$(NP_X) -DNP_Y=$(NP_Y) -DNP_T=$(NP_T)
	COMMS_LIB = comms.o
	FC = $(MPIFC)
else
	COMMS_LIB = uncomms.o
endif

ifeq ($(SITE_RANDOM), yes)
	RANDOM = site_random.o
	RANDOM_FLAGS = -DSITE_RANDOM
else
	RANDOM = random.o
	RANDOM_FLAGS = 
endif

default: bulk_rhmc compile_flags

.PHONY: clean

clean:
	rm -f bulk_rhmc *.o *.mod compile_flags

%.o: %.f90 Makefile
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} -c -o $@ $<

%.o: %.F90 Makefile
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} -c -o $@ $<

bulk_rhmc: params.o ${COMMS_LIB} ${RANDOM} bulk_rhmc_lib.o bulk_rhmc.f90
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} > $@
