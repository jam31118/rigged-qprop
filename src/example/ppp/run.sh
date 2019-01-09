#!/bin/bash
num_of_process=4
$QPROP_HOME/bin/imag-prop
$QPROP_HOME/bin/real-prop
$QPROP_DEP_DIR/openmpi/bin/mpiexec -n $num_of_process $QPROP_HOME/bin/ppp-mpi
$QPROP_DEP_DIR/openmpi/bin/mpiexec -n $num_of_process $QPROP_HOME/bin/eval-tsurff-mpi
bash $QPROP_HOME/src/util/data/combine_tsurff_data.sh

## NOTES
# The only diference between this example and the original examples
# is the presence of `use-ppp` parameter in `tsurff.param` file, 
# which configures `real-prop` program to stop propagation as soon as the laser field becomes off, 
# leaving the field-free post-progation to be done by `ppp`.

