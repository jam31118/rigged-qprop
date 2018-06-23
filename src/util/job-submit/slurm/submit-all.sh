#! /bin/bash

# Define submission scripts
script_submit_im_re="$QPROP_HOME/src/util/job-submit/slurm/imag-real.slurm"
script_submit_tsurff="$QPROP_HOME/src/util/job-submit/slurm/tsurff.slurm"

if [ ! -f $script_submit_im_re ]; then exit; fi
if [ ! -f $script_submit_tsurff ]; then exit; fi

# Submit imag and real propagation
mesg_im_re=$(sbatch $script_submit_im_re)
jid_im_re=$(echo $mesg_im_re | grep -oh "[1-9][0-9]*$")
mesg_tsurff=$(sbatch --dependency=afterok:$jid_im_re $script_submit_tsurff)
jid_tsurff=$(echo $mesg_tsurff | grep -oh "[1-9][0-9]*$")

