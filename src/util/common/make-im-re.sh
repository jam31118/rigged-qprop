#!/bin/bash

## Script log and error
#script_log_file=./make-im-re.log
#script_err_file=./make-im-re.err

### Initial logging
name_of_this_script="$0"
echo "[ LOG ] Name of this script : $name_of_this_script"

time_of_dispatch=$(date)
echo "[ LOG ] Dispatch Time of this script : $time_of_dispatch"

### Process input (don't need the number of process in make-im-re since there's no chance of mpirun)
#if [ -n "$1" ]; then num_of_process="$1"
#else (>&2 echo "[ERROR] No first argument \$1, the number of process."); exit 1; fi 

### Compile all the binaries
if [[ ! -f "Makefile" && ! -f "makefile" ]]; then (>&2 echo "No Makefile"); exit 1; fi
echo "[ LOG ] Starting compilation"
make all eval-tsurff-mpi > make.log 2>&1
exit_code="$?"
if [ "$exit_code" -ne 0 ]; then (>&2 echo "[ERROR] Failed to compile with exit code $exit_code"); exit 1; fi
echo "[ LOG ] Completed compilation"


### Basic Configuration


## Simulation program (binary) name list
im_bin=./hydrogen_im
re_bin=./hydrogen_re
#eval_bin=./eval-tsurff-mpi

## Parameter Files
initialParamFile="./initial.param"
#evalParamFile="./tsurff.param"

## Log files list
im_log="im.log"
re_log="re.log"
#eval_log="eval.log"


### Run simulation programs(binaries)
if [ ! -f "$im_bin" ]; then echo "[ERROR] No binary file $im_bin"; exit 1; fi
echo "[ LOG ] Started to run $im_bin"
$im_bin > $im_log 2>&1
exit_code="$?"
if [ "$exit_code" -ne 0 ]; then echo "[ERROR] $im_bin has exited with error code $exit_code, terminating the program"; exit 1; fi
echo "[ LOG ] $im_bin completed with time: $(date)"

if [ ! -f "$re_bin" ]; then echo "[ERROR] No binary file $re_bin"; exit 1; fi
echo "[ LOG ] Started to run $re_bin"
$re_bin > $re_log 2>&1
exit_code="$?"
if [ "$exit_code" -ne 0 ]; then echo "[ERROR] $re_bin has exited with error code $exit_code, terminating the program"; exit 1; fi
echo "[ LOG ] $re_bin completed with time: $(date)"



