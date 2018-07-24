#!/bin/bash

## Set Batch Options
#SBATCH --job-name=tsurff

#SBATCH --output=slurm.%x.%j.log

#SBATCH --ntasks=24
#SBATCH --nodes=2

#SBATCH --partition=g3,g4


## Run tsurff
hostname
$QPROP_DEP_DIR/openmpi/bin/mpiexec $QPROP_HOME/bin/eval-tsurff-mpi > eval.log 2>&1


## Combine partial data files produced by separate MPI processes
source $QPROP_HOME/src/util/data/combine_tsurff_data.sh

