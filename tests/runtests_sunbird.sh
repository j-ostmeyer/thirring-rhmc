#!/bin/bash

#SBATCH -t 5
#SBATCH --account=scw1000
#SBATCH -n 8
#SBATCH --qos=highpriority
#SBATCH --mem=20g


module purge
module load compiler/intel/2018/3
module load mpi/intel/2018/3

make runtests
