
# Do we need an extra rank for the MUST debugger?

include ./MkFlags
WITH_MUST=no
ifeq ($(WITH_MUST), yes)
	COMMS_FLAGS += -DWITH_MUST
endif


ifeq ($(MPI), yes)
	MPI_RUNNER = mpirun -n $$(($(NP_X) * $(NP_Y) * $(NP_T)))
else
	COMMS_LIB = uncomms.o
	MPI_RUNNER = 
endif



TESTS = test_dslash_1 test_dslash_3 test_dslash_5 test_dslashd_1 test_dslashd_3 \
       test_dslashd_5 test_dslash2d test_load test_save test_derivs test_qmrherm_0 \
       test_qmrherm_1 test_qmrherm_2 test_qmrherm_3 test_qmrherm_4 test_hamilton test_force \
       test_congrad test_congrad_1 test_measure test_gauss0 test_gaussp \
       test_halo_4 test_halo_4_real test_halo_5 test_halo_6 test_meson \
       test_dslashd_1_reqs test_dslashd_3_reqs test_dslashd_5_reqs \
       test_partitioning test_partitioning2 test_partitioning3 \
       test_comms_partitioning test_comms_partitioning2 test_dirac_split 

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
	rm -f $(TESTS) test_halo *.mod *.o compile_flags

.makefile.uptodate: Makefile ../MkRules ./MkFlags
	make clean
	touch .makefile.uptodate

params.o : test_params.F90 .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o params.o $<

test_utils.o : test_utils.F90 .makefile.uptodate comms.mod params.mod random.mod
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o test_utils.o $<


TOPDIR = ..

include ../MkRules

LIBOBJS += test_utils.o inverter_checks.o

MODS=$(patsubst %.o,%.mod,$(LIBOBJS))

.PRECIOUS : %.mod

inverter_checks.o  : inverter_checks.F90 \
        .makefile.uptodate dirac.mod trial.mod comms.mod
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o inverter_checks.o $<

test_%.o: test_%.F90 .makefile.uptodate test_utils.F90 test_utils.fh $(MODS)
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

test_%.o: test_%.f90 .makefile.uptodate test_utils.F90 test_utils.fh $(MODS)
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

$(TESTS) : % : %.o $(LIBOBJS) 
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -o $@ $^

