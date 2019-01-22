#!/bin/bash
#SBATCH --job-name lscpu
#SBATCH --output=slurm.%x.%j.log
#SBATCH --error=slurm.%x.%j.err

lscpu
