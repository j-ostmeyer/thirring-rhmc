
# HEAP OF OPTIONS THAT YOU MAY WANT TO TRY
#FC = /usr/lib64/mpi/gcc/mvapich2/bin/mpif90
#FC = /usr/lib64/mvapich2/bin/mpif90
#FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -heap-arrays -g
#FCFLAGS = -g -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1 #-heap-arrays -CB -traceback
#FCFLAGS = -g -O3 -march=native -mtune=native -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1
#FCFLAGS = -O0 -heap-arrays -warn all -C -traceback

COMPILER=gnu # either GCC or INTEL
MPI=no
NP_X=1
NP_Y=1
NP_T=1
SITE_RANDOM=no

#GNU SETTINGS
GNU_MPIFC    = mpif90
GNU_FC       = gfortran
GNU_FCFLAGS  = -O0 -Wall

#INTEL SETTINGS
INTEL_MPIFC  =mpiifort 
INTEL_FC     =ifort 
INTEL_FCFLAGS= -O0 -heap-arrays -warn all -C -traceback

ifeq ($(COMPILER), gnu)
	MPIFC  =$(GNU_MPIFC)
	FC     =$(GNU_FC)
	FCFLAGS=$(GNU_FCFLAGS)
else ifeq ($(COMPILER), INTEL)
	MPIFC  =$(INTEL_MPIFC)
	FC     =$(INTEL_FC)
	FCFLAGS=$(INTEL_FCFLAGS)
endif
$(info COMPILER: $(COMPILER))
$(info MPIFC   : $(MPIFC))
$(info FC      : $(FC))
$(info FCFLAGS : $(FCFLAGS))

ifeq ($(MPI), yes)
	COMMS_FLAGS = -DMPI -DNP_X=$(NP_X) -DNP_Y=$(NP_Y) -DNP_T=$(NP_T)
	COMMS_LIB = comms.o
	FC = $(MPIFC)
else ifeq ($(MPI), no)
	COMMS_LIB = uncomms.o
else 
    $(error MPI not correctly specified)
endif

ifeq ($(SITE_RANDOM), yes)
	RANDOM = site_random.o
	RANDOM_FLAGS = -DSITE_RANDOM
else ifeq ($(SITE_RANDOM), no)
	RANDOM = random.o
	RANDOM_FLAGS = 
else 
    $(error SITE_RANDOM not correctly specified)
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
