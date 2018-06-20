num_of_process=1
if [ -n "$1" ]; then num_of_process="$1"; fi
echo "num_of_process: $num_of_process" 

echo "[$(date)] imag-prop started"
$QPROP_HOME/bin/imag-prop > im.log 2>&1
echo "[$(date)] real-prop started"
$QPROP_HOME/bin/real-prop > re.log 2>&1
echo "[$(date)] eval-tsurff-mpi started"
$QPROP_DEP_DIR/openmpi/bin/mpiexec -np $num_of_process $QPROP_HOME/bin/eval-tsurff-mpi  > eval.log 2>&1

