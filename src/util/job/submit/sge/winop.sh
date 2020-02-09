# request "/bin/bash" as shell for job
#$ -S /bin/bash

# Import environment variable from where this job is submitted
#$ -V

# Use current directory as the working directory
#$ -cwd

# The Job name 
#$ -N winop

# Configure log files for stdout and stderr
#$ -o $JOB_NAME.$JOB_ID.log
#$ -e $JOB_NAME.$JOB_ID.err

#$ -pe mpich 1


date
echo "QPROP_HOME=${QPROP_HOME}"
echo "SGE_O_WORKDIR=${SGE_O_WORKDIR}"
echo "TMPDIR=${TMPDIR}"
echo "NSLOTS=${NSLOTS}"
echo "[ LOG ] machine list: "
cat $TMPDIR/machines


bin=winop

export LD_LIBRARY_PATH="$QPROP_DEP_DIR/openmpi/lib:$LD_LIBRARY_PATH"

echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

$QPROP_DEP_DIR/openmpi/bin/mpiexec --mca btl tcp,self -machinefile $TMPDIR/machines -np $NSLOTS $QPROP_HOME/bin/$bin > $bin.log 2>&1

echo "[ LOG ] ${bin} finished"

## Finalize data
#source $QPROP_HOME/src/util/data/combine_tsurff_data.sh
# [NOTE] needs to be fixed as MPI version
cat spectrum_polar0.dat > spectrum_polar.dat
cat spectrum_0.dat > spectrum.dat

echo "[ LOG ] ${bin} finished"

