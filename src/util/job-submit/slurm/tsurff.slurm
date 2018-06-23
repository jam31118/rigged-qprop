#!/bin/bash

## Set Batch Options
#SBATCH --job-name=tsurff

#SBATCH --output=slurm.%x.%j.log

#SBATCH --ntasks=12
#SBATCH --nodes=1

#SBATCH --partition=g1,g2,g3,g4


## Run tsurff
hostname
$QPROP_DEP_DIR/openmpi/bin/mpiexec $QPROP_HOME/bin/eval-tsurff-mpi


## Combine partial data files produced by separate MPI processes
source $QPROP_HOME/src/util/data/combine_tsurff_data.sh

