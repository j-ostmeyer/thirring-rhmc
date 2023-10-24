# Thirring Model in 2+1 dimensions with Domain Wall Fermions
## Simon Hands and Johann Ostmeyer

This tar-archive contains the program and data required to reproduce the results presented in *"Spectroscopy in the 2+1d Thirring Model with N=1 Domain Wall Fermions"*, [arXiv:2210.04790 [hep-lat]](https://arxiv.org/abs/2210.04790).

It contains a program-come-library for computations in lattice field theory with domain
wall fermions in 2+1D. Developed by Simon Hands, and refactored and
parallelised by Ed Bennett, Michele Mesiti, and the Swansea Academy of 
Advanced Computing RSE team. Features and analysis scripts added by Johann Ostmeyer.

The code is also available on github: [github.com/sa2c/thirring-rhmc](https://github.com/sa2c/thirring-rhmc).

## Data
The raw data produced by the program described below can be found in the directory `logruns`. It mirrors the structure on the DiRAC supercomputer.

The same data, but sorted in a systematic fashion, is repeated in the `correlators` directory. Therein it is analysed with the `ana.sh` and `ana_corr.R` scripts.

Data files:
- `fort.11` shows the average of the bose action (3rd column) for every trajectory (1st column).
- `fort.501` contains real and imaginary parts of the fermion correlator time series (time in 2nd column) in 3rd and 4th columns, defined by S_0 = 0.25*tr{gamma_0*S_f}  Cf. eqn (23).
- `fort.500` contains real and imaginary parts of the fermion correlator time series (time in 2nd column) in 3rd and 4th columns, defined by 0.25*tr{S_f} from which S_3 in eqn (23) is derived.
- `fort.302` contains meson correlator time series (time in 2nd column) involving gamma_5. Sum (difference) of 3rd and 4th columns is the gamma_5-Goldstone (gamma_5*gamma_3-Non-Goldstone).
- `fort.320` and `fort.321` contain meson correlator time series (time in 2nd column) involving the identity. 3rd and 4th columns are real and imaginary parts respectively. Sum (difference) of correlators in both files is the id-Goldstone (gamma_3-Non-Goldstone).
- `m_eff.csv` in the directory `correlators/results` contains all the fit results of the correlators using the appropriate effective mass ansatz.

## Analysis and Plot scripts
We use the `hadron` package in `R` for the analysis.

The bose action in `fort.11` is analysed for its autocorrelation time in order to extract an appropriate analysis range (avoid trajectories before thermalisation) and block size for the blocked bootstrap procedure. Then all the correlators are averaged, their errors estimated via bootstrap and effective masses calculated. The effective masses are fitted. These steps can be reproduced by running `./ana.sh -c` in the `correlators` directory. Intermediate results are stored as `.RData` files and individual correlators are plotted in the respective directories as `Rplots.pdf`. To summarise the results run `./ana.sh -p`.

The gnuplot script `plot_all.gp` automatically produces a summary of all the fit results. The resulting plots are located in `correlators/results` together with the script.

## Quick info / checklists for the main program
Remember to:
- checkout the `dev_mpithird` branch 
- Load the correct modules.

Regarding Modules, see `Machine-Specific options and caveats`.

Necessary files for:
- Compilation:
  - `MkFlags`: compiler choice, n of mpi ranks.
  - `params.F90`: `KSIZE`,`KSIZET`,`KTHIRD`
- Running:
  - compiled `bulk_rhmc` or link
  - `con` (optional, start from random if not present)
  - `midout`
  - `program_status` (optional)
  - `random_seed` (see `Random Seed settings`)
  - `remez` files:
    - `remez2`
    - `remez2g`
    - `remez4`
    - `remez4g`  

See examples.


## Compilation

Before compile time, the user may need to change the `MkFlags` file 
and the `params.F90` file. More details follow.

### Make-related files

Compilation is done using GNU `make`. In order to allow code reuse and 
modularity, the amount of information needed to compile is split in 3 files:
   * the `Makefile` itself: this is the file which is read by make, and there are 
different versions for it - a 'main' version, a 'test' version in the 'tests' 
directory, and a 'benchmark' version in the 'benchmarks' directory. This 
does not need to be modified (to my experience).
   * the `MkRules` file: the only version of this file is included in all makefiles
and contains the rules to build all the components of the program. 
This file does not need to be modified unless one wants to change 
compiler-specific flags, e.g. the flags for the compilers.
   * the `MkFlags` file: this is the file that the user is going to modify the most 
often. It contains the some details about the parallelization of the code. 
A possible content is the following:

```
COMPILER=INTEL# GNU,INTEL,CRAY,IBM,SPINTEL
MPI=yes#                                                                        
NP_X=2#                                                                         
NP_Y=2#                                                                         
NP_T=2#
NP_THIRD=2#
SITE_RANDOM=yes#
```

A brief explanation of the options follows.

#### Compiler choice
The choice of the compiler through the COMPILER variable just defines the 
compiler executable name and the flags to pass to the fortran compiler.
These are all defined in MkRules.

Notice that flags have been chosen carefully only for the intel compiler, 
and that there are no correct variables for the IBM compiler.


#### Compiling for parallel

Set `MPI` to `yes` to enable MPI. `NP_X`, `NP_Y`, `NP_T`, `NP_THIRD` must be set as the 
number of processes
in the x, y, t and third directions respectively; these are checked at compile time
so that the lattice dimensions are divisible by these, and at runtime against 
the number of ranks available. Also, at runtime, it is checked whether the local lattice
size along the third direction (`KTHIRD/NP_THIRD`) is divisible by 4, and if it isn't the execution aborts.
If `MPI` is set to `no`, `NP_X`, `NP_Y`, `NP_T` and `NP_THIRD` will be ignored.

#### Random numbers

Two options for random numbers are present: the original version, and a
site-based generator. The former guarantees the same random sequences as the
original code when running in serial, and works in parallel by re-seeding
other ranks based on the seed from rank 0 plus an offset - so results from
different parallelisations will use different random number sequences, so are
not guaranteed to be identical. The latter guarantees equivalence between runs
at different parallelisations, at the cost of losing equivalence with the
original version of the code even when running in serial.

##### Random Seed setting
The random seed can be set in a number of ways. It is read from/saved into the 
gauge configuration file (`iread` must be set to 1 for this to happen), 
but it can also be set manually in the code by also 
setting the `iseed` parameter in params.F90 to be different from zero.
A more flexible alternative is to write the seed in a file named `random_seed`. 
If that file is found, the new seed will be read from it and will replace the 
one read from the gauge configuration and the hard-coded one.
When the gauge configuration is saved, the seed is saved both in the 
gauge configuration file and in the 'random\_seed' file. If the user wants to 
change it, in order e.g. to have another statistically independent
 montecarlo history, it suffices to change the seed written in that file.

### params.F90 - Lattice volume and other parameters
The lattice extents, number of flavours, mass, control parameters for the
HMC (initial condition, whether to read or write, and how often), and the
number of iterations to use in the QMR inversion, are defined in `params.F90`.
This allows predefined numbers to be used in test cases, and removes the need
to edit the lines repeatedly throughout the code when e.g. the lattice volume 
changes.
Other version of the params.F90 file are used for tests and benchmarks - see 
the respective directories.

#### Caveats
* When modifying a Makefile or any file include in one, keep in mind that spaces,
even at the end of words, have a meaning. For example, setting a variable like this
```
MPI=yes # only characters after '#' will be ignored
```
is not the same as setting it like this
```
MPI=yes#
```

## Tests

Tests are found in the `tests` folder. They have their own `Makefile` that
works in exactly the same way as the one for the main code. Tests can be run
individually; no output means that they pass. Alternatively, `make runtests`
runs all tests.

The `.dat` files can be generated from the pre-MPI version of the code (or,
if you want, from the current version, but this should only be done when the
current version is known good), by setting the `generate` parameter in each program
to `.true.`.

## Benchmarks
In order to study the performance of the 'congrad' subroutine and of the various versions 
of the 'qmrherm' subroutine
(which make up most of the computational load of the program) two separate benchmark
programs have been created, with scripts and data needed to run them in many possible 
setups (e.g., benchsetup.py).
Moreover, a benchmark for the a full molecular dynamics trajectory has
also been written (it is not yet included in the auto-benchmark suite).
See the `benchmarks` directory.

## Machine-Specific options and caveats
* On the `enhpc` cluster it is necessary to load the `intel` and the `mvapich2`
modules. Moreover, the `INTEL_MPIFC` in the `MkRules` file must be set to 
`mpifort` instead of `mpiifort`.
* On the `hawk` and `sunbird` clusters the modules to load are `compiler/intel` and `mpi/intel`.
* On the `peta4` machine, the modules to load are 
  - `intel/bundles/complib/2018.4`
  It is better to unload the `2017.4` version 
  of the intel compiler and libraries
  first:
  ```
  module unload intel/bundles/complib/2017.4
  ```


