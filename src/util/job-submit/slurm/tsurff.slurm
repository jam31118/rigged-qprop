#!/bin/bash

#SBATCH --job-name=tsurff

#SBATCH --output=slurm.%x.%j.log
##SBATCH --error=error.%x.%j

#SBATCH --ntasks=20
##SBATCH --cpus-per-task=1
#SBATCH --nodes=2

##SBATCH --mem=2G

#SBATCH --partition=g1,g2,g3,g4

hostname
$QPROP_DEP_DIR/openmpi/bin/mpiexec $QPROP_HOME/bin/eval-tsurff-mpi
#$QPROP_DEP_DIR/openmpi/bin/mpiexec $QPROP_HOME/bin/eval-tsurff-mpi-intel

#srun bash $QPROP_HOME/src/util/run/tsurff.sh

