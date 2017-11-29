#!/bin/bash

# Set the name of the job
# (this gets displayed when you get a list of jobs on the cluster)
#SBATCH --job-name="sa2c_sjh_test"
#SBATCH -o sa2c_sjh_test.out
#SBATCH -e sa2c_sjh_test.err

# Specify the maximum wall clock time your job can use
# (Your job will be killed if it exceeds this)
#SBATCH --time=1:00:00

# Specify the amount of memory your job needs (in Mb)
# (Your job will be killed if it exceeds this for a significant length of time)
#SBATCH --mem-per-cpu=1024

# Specify the number of cpu cores your job requires
#SBATCH --ntasks=24

# Set up the environment
module load intel

# Run the application
export OMP_NUM_THREADS=$SLURM_NTASKS

if [ "$[OMP_NUM_THREADS]" -gt "1" ]
then
    FOLDER="`cat compile_flags` ${OMP_NUM_THREADS} threads"
else
    FOLDER="`cat compile_flags`"
fi

mkdir "$FOLDER"
cd "$FOLDER"
cp ../con ../midout ../remez* .
{ time ../bulk_rhmc; } 2> timings
rm con midout remez*
