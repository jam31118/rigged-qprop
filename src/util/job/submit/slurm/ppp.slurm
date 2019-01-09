#!/bin/bash

## Set Batch Options
#SBATCH --job-name=ppp

#SBATCH --output=slurm.%x.%j.log

#SBATCH --ntasks=12
#SBATCH --nodes=1

#SBATCH --partition=g3,g4


## Run tsurff
hostname
bin="ppp-mpi"
log_file="${bin}.log"
echo "log file name: $log_file"
echo "binary file name: $bin"
$QPROP_DEP_DIR/openmpi/bin/mpiexec $QPROP_HOME/bin/$bin > $log_file 2>&1

