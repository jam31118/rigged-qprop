#!/bin/bash

# request 'bash' as shell for job
#$ -S /bin/bash

# Import environment variable from where this job is submitted
#$ -V

# The Job name 
#$ -N single

# Configure log files for stdout and stderr
#$ -o $JOB_NAME.$JOB_ID.log
#$ -e $JOB_NAME.$JOB_ID.err

cd $SGE_O_WORKDIR

log_file="${JOB_NAME}.${JOB_ID}.log"
date > $log_file

bin=imag-prop
$QPROP_HOME/bin/$bin > $bin.log 2>&1
echo "[ LOG ] ${bin} finished" >> $log_file

bin=real-prop
$QPROP_HOME/bin/$bin > $bin.log 2>&1
echo "[ LOG ] ${bin} finished" >> $log_file

