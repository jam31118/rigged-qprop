#!/bin/bash

#SBATCH --job-name=imag-real

#SBATCH --output=slurm.%x.%j.log

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --partition=g3,g4

bash $QPROP_HOME/src/util/run/imag.sh
bash $QPROP_HOME/src/util/run/real.sh

