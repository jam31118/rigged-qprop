#! /bin/bash

## Source required scripts
source $QPROP_HOME/src/util/shell/colors.sh

## Get sbatchargument for tsurff calculation
tsurff_submit_arg="$@"
echo -e "${LOG} Given tsurff_submit_arg: $tsurff_submit_arg"

## Define submission scripts
SCRIPT_DIR="$QPROP_HOME/src/util/job/submit/slurm"
script_submit_im_re="$SCRIPT_DIR/imag-real.slurm"
script_submit_tsurff="$SCRIPT_DIR/tsurff.slurm"

## Check whether the script files exist
if [ ! -f $script_submit_im_re ]; then (2>&1 echo -e "${ERROR} Could not find script: $script_submit_im_re"); exit 1; fi
if [ ! -f $script_submit_tsurff ]; then (2>&1 echo -e "${ERROR} Could not fine script: $script_submit_tsurff"); exit 1; fi

## Submit imag and real propagation
# submit imag-real
mesg_im_re=$(sbatch $script_submit_im_re)
if [ "$?" -ne 0 ]; then (2>&1 echo -e "${ERROR} Failed to submit $script_submit_im_re"); fi
jid_im_re=$(echo $mesg_im_re | grep -oh "[1-9][0-9]*$")
echo -e "${OK} [$(date)] Submitted $script_submit_im_re with job id: $jid_im_re"

# submit tsurff
mesg_tsurff=$(sbatch $tsurff_submit_arg --dependency=afterok:$jid_im_re $script_submit_tsurff)
if [ "$?" -ne 0 ]; then (2>&1 echo -e "${ERROR} Failed to submit $script_submit_tsurff"); fi
jid_tsurff=$(echo $mesg_tsurff | grep -oh "[1-9][0-9]*$")
echo -e "${OK} [$(date)] Submitted $script_submit_tsurff with job id: $jid_tsurff"

