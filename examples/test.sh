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


for file in con random_seed program_status
do 
  cp $file $file.$SLURM_JOB_ID
done


BOOKKEEPING=bookkeeping.txt

if [ ! -f "$BOOKKEEPING" ]
then
    echo '# Line count' >> $BOOKKEEPING
    echo UnixTimeStamp MeasType SLURM_JOB_ID fort.100 fort.11 fort.200 control output>> $BOOKKEEPING
fi

echo $(date +"%s") start $SLURM_JOB_ID $(cat fort.100  2>/dev/null| wc -l ) \
 $(cat fort.11 2>/dev/null | wc -l) $( cat fort.200 2>/dev/null | wc -l ) \
 $(cat control 2>/dev/null | wc -l) $(cat output 2>/dev/null | wc -l) >> $BOOKKEEPING

echo "Slurm Job ID: $SLURM_JOB_ID" >> output
$MPIRUNNER  ./bulk_rhmc

echo $(date +"%s") end $SLURM_JOB_ID $(cat fort.100  2>/dev/null| wc -l ) \
 $(cat fort.11 2>/dev/null | wc -l) $( cat fort.200 2>/dev/null | wc -l ) \
 $(cat control 2>/dev/null | wc -l) $(cat output 2>/dev/null | wc -l) >> $BOOKKEEPING


