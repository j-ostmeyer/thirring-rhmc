# dwf3d_sjh / bulk_rhmc

A program-come-library for computations in lattice field theory with domain
wall fermions in 2+1D. Developed by Simon Hands, and refactored and
parallelised by Ed Bennett and the Swansea Academy of Advanced Computing RSE
team.

## Compilation

Compilation is done using GNU `make`. In the simplest case, set a Fortran
compiler as `FC` and the flags to pass to it as `FCFLAGS`.

### Compiling in parallel

Set `MPIFC` to the MPI wrapped Fortran compiler, and set `MPI` to `yes` to
enable MPI. `NP_X`, `NP_Y`, and `NP_T` must be set as the number of processes
in the x, y, and t directions respectively; these are checked at compile time
so that the lattice dimensions are divisible by these, and at runtime against 
the number of ranks available. 

### Random numbers

Two options for random numbers are present: the original version, and a
site-based generator. The former guarantees the same random sequences as the
original code when running in serial, and works in parallel by re-seeding
other ranks based on the seed from rank 0 plus an offset - so results from
different parallelisations will use different random number sequences, so are
not guaranteed to be identical. The latter guarantees equivalence between runs
at different parallelisations, at the cost of losing compatibility with the
original version of the code even when running in serial.

### Lattice volume and other parameters

The lattice extents, number of flavours, mass, control parameters for the
HMC (initial condition, whether to read or write, and how often), and the
number of iterations to use in the QMR inversion, are defined in `params.F90`.
This allows predefined numbers to be used in test cases, and removes the need
to edit the lines repeatedly throughout the code when e.g. the lattice volume 
changes.

## Tests

Tests are found in the `tests` folder. They have their own `Makefile` that
works in exactly the same way as the one for the main code. Tests can be run
individually; no output means that they pass. Alternatively, `make runtests`
runs all tests.
