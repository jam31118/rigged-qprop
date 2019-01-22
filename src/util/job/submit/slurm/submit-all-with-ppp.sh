#! /bin/bash


## Source required scripts
source $QPROP_HOME/src/util/shell/colors.sh


## Get sbatchargument for tsurff calculation
parallel_submit_arg="$@"
echo -e "${LOG} Given parallel_submit_arg: $parallel_submit_arg"


## Define submission scripts
SCRIPT_DIR="$QPROP_HOME/src/util/job/submit/slurm"
script_submit_im_re="$SCRIPT_DIR/imag-real.slurm"
script_submit_tsurff="$SCRIPT_DIR/tsurff.slurm"
script_submit_ppp_mpi="$SCRIPT_DIR/ppp.slurm"


## Check whether the script files exist
for script in $script_submit_im_re $script_submit_tsurff $script_submit_ppp_mpi
do
  if [ ! -f $script ]; then (2>&1 echo -e "${ERROR} Could not find script: $script"); exit -1; fi
done
#if [ ! -f $script_submit_im_re ]; then (2>&1 echo -e "${ERROR} Could not find script: $script_submit_im_re"); exit 1; fi
#if [ ! -f $script_submit_tsurff ]; then (2>&1 echo -e "${ERROR} Could not fine script: $script_submit_tsurff"); exit 1; fi


## Submit imag and real propagation
# submit imag-real
mesg_im_re=$(sbatch $script_submit_im_re)
if [ "$?" -ne 0 ]; then (2>&1 echo -e "${ERROR} Failed to submit $script_submit_im_re"); fi
jid_im_re=$(echo $mesg_im_re | grep -oh "[1-9][0-9]*$")
echo -e "${OK} [$(date)] Submitted $script_submit_im_re with job id: $jid_im_re"

# submit ppp-mpi
mesg_ppp_mpi=$(sbatch $parallel_submit_arg --dependency=afterok:$jid_im_re $script_submit_ppp_mpi)
if [ "$?" -ne 0 ]; then (2>&1 echo -e "${ERROR} Failed to submit $script_submit_mpi"); fi
jid_ppp_mpi=$(echo $mesg_ppp_mpi | grep -oh "[1-9][0-9]*$")
echo -e "${OK} [$(date)] Submitted $script_submit_ppp_mpi with job id: $jid_ppp_mpi"

# submit tsurff
mesg_tsurff=$(sbatch $parallel_submit_arg --dependency=afterok:$jid_ppp_mpi $script_submit_tsurff)
if [ "$?" -ne 0 ]; then (2>&1 echo -e "${ERROR} Failed to submit $script_submit_tsurff"); fi
jid_tsurff=$(echo $mesg_tsurff | grep -oh "[1-9][0-9]*$")
echo -e "${OK} [$(date)] Submitted $script_submit_tsurff with job id: $jid_tsurff"

