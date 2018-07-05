
# HEAP OF OPTIONS THAT YOU MAY WANT TO TRY
#FC = /usr/lib64/mpi/gcc/mvapich2/bin/mpif90
#FC = /usr/lib64/mvapich2/bin/mpif90
#FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -heap-arrays -g
#FCFLAGS = -g -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1 #-heap-arrays -CB -traceback
#FCFLAGS = -g -O3 -march=native -mtune=native -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1
#FCFLAGS = -O0 -heap-arrays -warn all -C -traceback

COMPILER=INTEL# either GNU or INTEL
MPI=yes
NP_X=2
NP_Y=2
NP_T=2
SITE_RANDOM=yes

#GNU SETTINGS
GNU_MPIFC    = mpif90
GNU_FC       = gfortran
GNU_FCFLAGS  = -O0 -Wall -ffree-line-length-none

#INTEL SETTINGS
INTEL_MPIFC  =mpiifort 
INTEL_FC     =ifort 
INTEL_FCFLAGS= -O3 -heap-arrays

ifeq ($(COMPILER), GNU)
MPIFC  =$(GNU_MPIFC)
FC     =$(GNU_FC)
FCFLAGS=$(GNU_FCFLAGS)
else ifeq ($(COMPILER), INTEL)
MPIFC  =$(INTEL_MPIFC)
FC     =$(INTEL_FC)
FCFLAGS=$(INTEL_FCFLAGS)
else 
$(error COMPILER not correctly specified (watch for whitespaces)) 
endif
$(info COMPILER: $(COMPILER))
$(info MPIFC   : $(MPIFC))
$(info FC      : $(FC))
$(info FCFLAGS : $(FCFLAGS))

OBJS = bulk_rhmc.o avgitercounts.o dirac.o dum1.o dwf3d_lib.o \
       gauge.o gaussian.o gforce.o params.o phizero.o \
       qmrherm_scratch.o remez.o remezg.o trial.o vector.o\
       comms.o random.o

default: bulk_rhmc compile_flags

.PHONY: clean

COMPILE = $(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} -c 



bulk_rhmc.o : bulk_rhmc.f90 Makefile dwf3d_lib.mod
	$(COMPILE) -o $*.o $<

avgitercounts.o avgitercounts.mod : avgitercounts.F90 Makefile
	echo $(COMPILE)
	$(COMPILE) -o $*.o $<

dirac.o dirac.mod : dirac.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

dum1.o dum1.mod : dum1.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

dwf3d_lib.o dwf3d_lib.mod : dwf3d_lib.F90 Makefile avgitercounts.mod \
    comms.mod dirac.mod dum1.mod gauge.mod gaussian.mod gforce.mod \
    params.mod phizero.mod qmrherm_scratch.mod random.mod remez.mod \
    remezg.mod trial.mod vector.mod
	$(COMPILE) -o $*.o $<

gauge.o gauge.mod : gauge.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

gaussian.o gaussian.mod : gaussian.F90 Makefile comms.mod params.mod random.mod 
	$(COMPILE) -o $*.o $<

gforce.o gforce.mod : gforce.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

params.o params.mod : params.F90 Makefile
	$(COMPILE) -o $*.o $<

phizero.o phizero.mod   : phizero.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

qmrherm_scratch.o qmrherm_scratch.mod : qmrherm_scratch.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

remez.o remez.mod : remez.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

remezg.o remezg.mod : remezg.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

trial.o trial.mod : trial.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

vector.o vector.mod : vector.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<


# COMMUNICATION-RELATED 
ifeq ($(MPI), yes)
COMMS_FLAGS = -DMPI -DNP_X=$(NP_X) -DNP_Y=$(NP_Y) -DNP_T=$(NP_T)
FC = $(MPIFC)

comms.o comms.mod: comms.F90 Makefile params.mod
	$(COMPILE) -o $*.o $<

else ifeq ($(MPI), no)

comms.o comms.mod: uncomms.f90 Makefile params.mod
	$(COMPILE) -o $*.o $<

else 
$(error MPI not correctly specified (watch for whitespaces))
endif

# RNG-RELATED
ifeq ($(SITE_RANDOM), yes)
RANDOM_FLAGS = -DSITE_RANDOM

random.o random.mod: site_random.f90 Makefile comms.mod params.mod
	$(COMPILE) -o $*.o $<

else ifeq ($(SITE_RANDOM), no)

RANDOM_FLAGS = 

random.o random.mod: random.f90 Makefile comms.mod
	$(COMPILE) -o $*.o $<
else 
    $(error SITE_RANDOM not correctly specified (watch for whitespaces))
endif


clean:
	rm -f bulk_rhmc *.o *.mod compile_flags


bulk_rhmc: $(OBJS)
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} > $@

