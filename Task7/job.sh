#!/bin/sh
#module add mpi/openmpi-local
module add openmpi
for ((i = 0; i < 7; i++))
do
echo $i
mpirun -np $((2**i)) --oversubscribe ./main
done
