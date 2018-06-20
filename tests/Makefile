
# HEAP OF OPTIONS THAT YOU MAY WANT TO TRY
#FC = /usr/lib64/mpi/gcc/mvapich2/bin/mpif90B
#FC = /usr/lib64/mvapich2/bin/mpif90
#MPIFC = mpiifort
#FC = ifort
#FCFLAGS = -no-prec-div -fp-model fast=2 -xHost -O3 -check_mpi -g
#FCFLAGS = -g -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1 #-heap-arrays -CB -traceback
#FCFLAGS = -g -O3 -march=native -mtune=native -DMPI -DNP_X=1 -DNP_Y=1 -DNP_T=1
#FCFLAGS = -g -O0 -CB -warn all #-traceback

COMPILER=INTEL# either GNU or INTEL
MPI=no
NP_X=1
NP_Y=1
NP_T=1
SITE_RANDOM=yes
# Do we need an extra rank for the MUST debugger?
WITH_MUST=no


#GNU SETTINGS
GNU_MPIFC    = mpif90
GNU_FC       = gfortran
GNU_FCFLAGS  = -O0 -Wall -ffree-line-length-none

#INTEL SETTINGS
INTEL_MPIFC  =mpiifort 
INTEL_FC     =ifort 
INTEL_FCFLAGS= -O0 -g -heap-arrays -warn all -C -traceback

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


ifeq ($(MPI), yes)
	COMMS_FLAGS = -DMPI -DNP_X=$(NP_X) -DNP_Y=$(NP_Y) -DNP_T=$(NP_T)
	COMMS_LIB = comms.o
	MPI_RUNNER = mpirun -n $$(($(NP_X) * $(NP_Y) * $(NP_T)))
	FC = $(MPIFC)
else
	COMMS_LIB = uncomms.o
	MPI_RUNNER = 
endif

ifeq ($(SITE_RANDOM), yes)
	RANDOM = site_random.o
	RANDOM_FLAGS = -DSITE_RANDOM
else
	RANDOM = random.o
	RANDOM_FLAGS = 
endif

ifeq ($(WITH_MUST), yes)
	COMMS_FLAGS += -DWITH_MUST
endif

TESTS = test_dslash_1 test_dslash_3 test_dslash_5 test_dslashd_1 test_dslashd_3 \
       test_dslashd_5 test_dslash2d test_load test_save test_derivs test_qmrherm_0 \
       test_qmrherm_1 test_qmrherm_2 test_qmrherm_3 test_qmrherm_4 test_hamilton test_force \
       test_congrad test_congrad_1 test_measure test_gauss0 test_gaussp \
       test_halo_4 test_halo_4_real test_halo_5 test_halo_6 test_meson

default: $(TESTS)
tests: $(TESTS)

runtests: $(TESTS)
	@for test in $(TESTS); do \
	    printf "$$test: "; \
	    OUTPUT=$$($(MPI_RUNNER) ./$$test); \
	    if [ -z "$$OUTPUT" ]; then \
	        echo "OK"; \
	    else \
	        echo; \
	        echo $$OUTPUT; \
	    fi; \
	done

clean:
	rm -f $(TESTS) test_halo *.mod *.o

.makefile.uptodate: Makefile
	make clean
	touch .makefile.uptodate

%.o: ../%.f90 test_params.o .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

%.o: ../%.F90 test_params.o .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

bulk_rhmc_lib.o : ../bulk_rhmc_lib.F90 test_params.o $(COMMS_LIB) $(RANDOM) .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

site_random.o random.o : %.o : ../%.f90 test_params.o $(COMMS_LIB) .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

inverter_checks.o : inverter_checks.F90  .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<


test_%.o: test_%.f90 .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

test_params.o: test_params.F90 .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

test_utils.o: test_utils.F90 $(COMMS) $(RANDOM) .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

test_halo_%: test_params.o $(COMMS_LIB) test_halo_%.F90
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -o $@ $^

$(TESTS) : % : test_params.o $(COMMS_LIB) $(RANDOM) bulk_rhmc_lib.o \
          test_utils.o inverter_checks.o %.F90
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -o $@ $^

compile_flags:
	echo $(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) > $@
