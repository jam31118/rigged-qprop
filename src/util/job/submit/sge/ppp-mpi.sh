#!/bin/bash

# request 'bash' as shell for job
#$ -S /bin/bash

# Import environment variable from where this job is submitted
#$ -V

# Use current directory as the working directory
#$ -cwd

# The Job name 
#$ -N ppp-mpi

# Configure log files for stdout and stderr
#$ -o $JOB_NAME.$JOB_ID.log
#$ -e $JOB_NAME.$JOB_ID.err

#$ -pe mpich 1


log_file="${JOB_NAME}.${JOB_ID}.log"
date > $log_file

export LD_LIBRARY_PATH="$QPROP_DEP_DIR/openmpi/lib:$LD_LIBRARY_PATH"

bin=ppp-mpi
$QPROP_DEP_DIR/openmpi/bin/mpiexec --mca btl tcp,self -machinefile $TMPDIR/machines -np $NSLOTS $QPROP_HOME/bin/$bin > $bin.log 2>&1
echo "[ LOG ] ${bin} finished" >> $log_file

