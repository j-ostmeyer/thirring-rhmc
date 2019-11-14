#!/bin/bash

# This is a simple script that will find and run all the benchmarks 
# that were created by the script benchsetup.py. Refer to the docs in it
# for details.

# Settings may be SUNBIRD specific.

#SBATCH -N 1
#SBATCH --partition=compute
#SBATCH --time=20
#SBATCH --account=scw1000
#SBATCH --exclusive

module load compiler/intel/2018/3
module load mpi/intel/2018/3
for script in $(find  benchmarks* -name scriptrun1)
do
   echo "found /$script"
   DIR=$(dirname $script)
   EXE=$(basename $script)
   echo cd $DIR
   cd $DIR
   echo bash ./$EXE
   bash ./$EXE
   cd -
done
