#!/bin/bash

#SBATCH -t 60
#SBATCH --account=scw1000
##SBATCH -n 40
#SBATCH -o test%J
#SBATCH --qos=highpriority
#SBATCH --mem=20g
#SBATCH --oversubscribe
#SBATCH --job-name="BENCHMARK"
#SBATCH --dependency=singleton

module load mpi/intel/2018/3 compiler/intel/2018/3
MPIRUNNER='mpirun -n 8' 

echo "Slurm Job ID: $SLURM_JOB_ID" >> output

cp con con.$SLURM_JOB_ID
cp random_seed random_seed.$SLURM_JOB_ID


BOOKKEEPING=bookkeeping.txt

if [ ! -f "$BOOKKEEPING" ]
then
    echo '# Line count' >> $BOOKKEEPING
    echo SLURM_JOB_ID fort.100 fort.11 fort.200 control >> $BOOKKEEPING
fi

echo $SLURM_JOB_ID $(cat fort.100 | wc -l ) $(cat fort.11 | wc -l) \
    $( cat fort.200 | wc -l ) $(cat control | wc -l) >> $BOOKKEEPING


$MPIRUNNER  ./bulk_rhmc
