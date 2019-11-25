#!/bin/bash

num_slot=1
if [ -n "$1" ]; then num_slot="$1"; fi

echo "number of slot(s): $num_slot"

jid_imag_real=$(qsub -terse $QPROP_HOME/src/util/job/submit/sge/single.sh)
echo jid_imag_real: $jid_imag_real
jid_ppp_mpi=$(qsub -terse -pe mpich $num_slot -hold_jid $jid_imag_real $QPROP_HOME/src/util/job/submit/sge/ppp-mpi.sh)
echo jid_ppp_mpi: $jid_ppp_mpi
jid_tsurff_mpi=$(qsub -terse -pe mpich $num_slot -hold_jid $jid_imag_real,$jid_ppp_mpi $QPROP_HOME/src/util/job/submit/sge/tsurff.sh)
echo jid_tsurff_mpi: $jid_tsurff_mpi

