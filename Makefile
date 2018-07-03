
# HEAP OF OPTIONS THAT YOU MAY WANT TO TRY
#FC = /usr/lib64/mpi/gcc/mvapich2/bin/mpif90
#FC = /usr/lib64/mvapich2/bin/mpif90
#FCFLAGS = -ipo -no-prec-div -fp-model fast=2 -xHost -O3 -heap-arrays -g
#FCFLAGS = -g -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1 #-heap-arrays -CB -traceback
#FCFLAGS = -g -O3 -march=native -mtune=native -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1
#FCFLAGS = -O0 -heap-arrays -warn all -C -traceback

COMPILER=INTEL# either GNU or INTEL
MPI=no
NP_X=1
NP_Y=1
NP_T=1
SITE_RANDOM=yes

#GNU SETTINGS
GNU_MPIFC    = mpif90
GNU_FC       = gfortran
GNU_FCFLAGS  = -O0 -Wall -ffree-line-length-none

#INTEL SETTINGS
INTEL_MPIFC  =mpiifort 
INTEL_FC     =ifort 
INTEL_FCFLAGS= -O0 -g -heap-arrays -warn all

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


OBJS = bulk_rhmc.o avgitercounts.o dirac.o dum1.o dwf3d_lib.o gauge.o gforce.o params.o phizero.o qmrherm_only.o qmrherm_scratch.o remez.o remezg.o trial.o vector.o comms.o random.o

DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
ifeq ($(COMPILER), GNU)
	DEPFLAGS = -cpp -MMD -E -MP -MF $(DEPDIR)/$*.Td
else ifeq ($(COMPILER), INTEL)
	DEPFLAGS = -module -P -syntax-only -gen-dep=$(DEPDIR)/$*.Td
else 
    $(error COMPILER not correctly specified (watch for whitespaces))
endif


default: bulk_rhmc compile_flags

.PHONY: clean

COMPILE = $(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} -c 
POSTGETDEPS = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d 
GETDEPS = $(FC) $(DEPFLAGS) $(FCFLAGS)  $(COMMS_FLAGS) ${RANDOM_FLAGS} -c 

%.d : %.f90 Makefile
	$(GETDEPS) $<
	$(POSTGETDEPS) 

%.o : %.f90 $(DEPDIR)/%.d Makefile
	$(COMPILE) -o $*.o $<

%.d : %.F90 Makefile
	$(GETDEPS) $<
	$(POSTGETDEPS) 

%.o : %.F90 $(DEPDIR)/%.d Makefile
	$(COMPILE) -o $*.o $<


ifeq ($(MPI), yes)
	COMMS_FLAGS = -DMPI -DNP_X=$(NP_X) -DNP_Y=$(NP_Y) -DNP_T=$(NP_T)
	FC = $(MPIFC)
comms.d : comms.F90 Makefile
	$(GETDEPS) $<
	$(POSTGETDEPS) 
        
comms.o comms.mod: comms.F90 $(DEPDIR)/comms.d Makefile
	$(COMPILE) -o $*.o $<
else ifeq ($(MPI), no)

comms.o comms.mod: comms.F90 Makefile
	$(GETDEPS) $<
	$(POSTGETDEPS) 
comms.o comms.mod: uncomms.f90 $(DEPDIR)/comms.d Makefile
	$(COMPILE) -o $*.o $<
else 
    $(error MPI not correctly specified (watch for whitespaces))
endif

ifeq ($(SITE_RANDOM), yes)
	RANDOM_FLAGS = -DSITE_RANDOM
random.o random.mod: site_random.f90 Makefile
random.o random.mod: site_random.f90 $(DEPDIR)/random.d Makefile
	$(COMPILE) -o $*.o $<
else ifeq ($(SITE_RANDOM), no)
	RANDOM_FLAGS = 
random.o random.mod: random.f90 Makefile
random.o random.mod: random.f90 $(DEPDIR)/random.d Makefile
	$(COMPILE) -o $*.o $<
else 
    $(error SITE_RANDOM not correctly specified (watch for whitespaces))
endif


$(DEPDIR)/%.d: ; # generated rules are an actual make target
.PRECIOUS: $(DEPDIR)/%.d

# including generated rules 
include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(OBJS))))
	
clean:
	rm -f bulk_rhmc *.o *.mod compile_flags && rm -r .d


bulk_rhmc: $(OBJS)
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} > $@


