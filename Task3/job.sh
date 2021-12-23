#!/bin/sh
#module add mpi/openmpi-local
module add openmpi
mpirun -np 2 ./main
