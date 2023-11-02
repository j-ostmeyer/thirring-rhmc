#!/bin/bash -l

# Batch script to run an MPI parallel job under SGE with Intel MPI.

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=2:00:0

# Request 1 gigabyte of RAM per process (must be an integer followed by M, G, or T)
#$ -l mem=2G

# Set the name of the job.
#$ -N ThirringTest

# Select the MPI parallel environment and 16 processes.
#$ -pe mpi 32

# Set the working directory to somewhere in your scratch space.
# Replace "<your_UCL_id>" with your UCL user ID :
#$ -cwd

# Run our MPI job.  GERun is a wrapper that launches MPI jobs on our clusters.
date
gerun /home/ccaecai/Scratch/thirring-rhmc/tests/e2e_tests/<EXECUTABLE_PATH>
date
