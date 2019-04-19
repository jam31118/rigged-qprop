#!/bin/bash

# request "/bin/sh" as shell for job
#$ -S /bin/sh

# Import environment variable from where this job is submitted
#$ -V

# The Job name 
#$ -N single

#$ -o $HOME/job-log/$JOB_NAME.$JOB_ID.log
#$ -e $HOME/job-log/$JOB_NAME.$JOB_ID.err

cd $SGE_O_WORKDIR

#date > log.log
#echo "QPROP_HOME=${QPROP_HOME}" >> log.log
#echo "SGE_O_WORKDIR=${SGE_O_WORKDIR}"

log_file="${JOB_NAME}.${JOB_ID}.log"
date > $log_file

bin=imag-prop
$QPROP_HOME/bin/$bin > $bin.log 2>&1
echo "[ LOG ] ${bin} finished" >> $log_file

bin=real-prop
$QPROP_HOME/bin/$bin > $bin.log 2>&1
echo "[ LOG ] ${bin} finished" >> $log_file

# SGE_O_WORKDIR/exec_file
