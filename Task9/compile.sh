#!/bin/sh
#module add mpi/openmpi-local
module add openmpi
mpic++ main.cpp -o main
mpic++ main1.cpp -o main1
