#!/bin/bash

#SBATCH --job-name=imag-real

#SBATCH --output=slurm.%x.%j.log

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --partition=g1,g2,g3,g4

srun bash $QPROP_HOME/src/util/run/imag.sh
srun bash $QPROP_HOME/src/util/run/real.sh

