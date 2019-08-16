
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



BENCHMARKS = benchmark_congrad benchmark_qmrherm_1 \
 benchmark_qmrherm_1_sp benchmark_qmrherm_split1 \
 benchmark_qmrherm_split_nodir1 benchmark_full_md

default: $(BENCHMARKS)
benchmarks: $(BENCHMARKS)

runbenchmarks: $(BENCHMARKS)
	@for benchmark in $(BENCHMARKS); do \
	    printf "$$benchmark: "; \
	    OUTPUT=$$($(MPI_RUNNER) ./$$benchmark); \
	    if [ -z "$$OUTPUT" ]; then \
	        echo "OK"; \
	    else \
	        echo; \
	        echo $$OUTPUT; \
	    fi; \
	done


clean:
	rm -f $(BENCHMARKS) benchmark_halo *.mod *.o compile_flags

.makefile.uptodate: Makefile ../MkRules ./MkFlags
	make clean
	touch .makefile.uptodate

params.o : benchmark_params.F90 .makefile.uptodate
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o params.o $<

TOPDIR = ..

include ../MkRules

MODS=$(patsubst %.o,%.mod,$(LIBOBJS))

.PRECIOUS : %.mod

benchmark_%.o: benchmark_%.F90 .makefile.uptodate $(MODS)
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -c -o $@ $<

$(BENCHMARKS) : % : %.o $(LIBOBJS) 
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) $(RANDOM_FLAGS) -dynamic -o $@ $^

