#!/bin/bash

source $QPROP_HOME/src/util/shell/colors.sh


# Process command-line arguments

## for selecting binary file
bin_name=ppp
if [ -n "$1" ]; then bin_name="$1"; fi
echo "${LOG} binary file name: '$bin_name'"

## for setting num_of_process
in_parallel=false
num_of_process=1
if [ -n "$2" ]; then num_of_process="$2"; if [ "$bin_name" == "ppp-mpi" ]; then in_parallel=true; fi; fi



# Compilation
make -j4 $bin_name -C $QPROP_HOME/src/main
make_success="$?"
if [ "$make_success" -ne "0" ]; then (2>&1 echo "${ERROR} during 'make'"); exit -1; fi

## Copy field-present propagation data to the post-propagation directroy
## .. if there is already one, erase it first.
if [ -d "ati-ppp" ]; then rm -rf ati-ppp; fi
cp -r ati-before-ppp ati-ppp

## Enter the calculation directory
cd ati-ppp

## Run calculation
if $in_parallel
then $QPROP_DEP_DIR/openmpi/bin/mpiexec -n $num_of_process $QPROP_HOME/src/main/$bin_name; ppp_success="$?"
else $QPROP_HOME/src/main/$bin_name; ppp_success="$?"
fi

## Go back to the original directory
cd ..
if [ "$ppp_success" -ne "0" ]; then (2>&1 echo "${ERROR} failed to run '$bin_name' successfully"); exit -1; fi



# Run comparison scripts
source ~/py/venv/sci/bin/activate
echo "${LOG} running comparison script for wavefunction files"
python ./compare-wf.py
echo "${LOG} running comparison script for *.raw files"
python ./compare-raw.py

