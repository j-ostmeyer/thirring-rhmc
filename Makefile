
include ./MkFlags
.makefile.uptodate: Makefile MkRules MkFlags
	make clean
	touch .makefile.uptodate

default: bulk_rhmc compile_flags

# PARAMETERS
params.o params.mod : params.F90 .makefile.uptodate
	$(COMPILE) -o params.o params.F90

TOPDIR = .
include ./MkRules

.PRECIOUS : %.mod
# MAIN PROGRAM
bulk_rhmc: $(LIBOBJS) bulk_rhmc.o
	$(FC) $(FCFLAGS) $(COMMS_FLAGS) ${RANDOM_FLAGS} -o $@ $^

clean:
	rm -f bulk_rhmc *.o *.mod compile_flags *_sp.F90


