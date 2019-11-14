#!/bin/bash
TOTP=5
for i in `seq $TOTP` 
do
(mpirun -n $((40/TOTP)) ./benchmark_qmrherm_1 &> out$i)&
echo 'Launched ' $i
done 
wait
