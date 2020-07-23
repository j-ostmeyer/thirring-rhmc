#!/bin/bash
#SBATCH --job-name=benchmark
#SBATCH --time=20
#SBATCH --mem-per-cpu=1024
#SBATCH --exclusive=user
#SBATCH --ntasks=SEDNPROCS
#SBATCH --account=DIRAC-DP137-SL4-CPU
#SBATCH --partition=skylake

# Set up the environment
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load slurm rhel7/global rhel7/default-peta4 intel/bundles/complib/2018.4

MPIRUNNER='mpirun -n SEDNPROCS' 

time $MPIRUNNER  ./benchmark_congrad
time $MPIRUNNER  ./benchmark_qmrherm_1
time $MPIRUNNER  ./benchmark_qmrherm_1_sp

