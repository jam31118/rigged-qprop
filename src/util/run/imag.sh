#!/bin/bash

source $QPROP_HOME/src/util/shell/colors.sh

bin_name="imag-prop"
bin_path=$QPROP_HOME/bin/$bin_name
log_name="im.log"

if [ ! -f "$bin_path" ]; then (2>&1 echo "${ERROR} Cannot find $bin_path"); exit -1; fi 
echo "${LOG}[$(date)] $bin_name started"
$bin_path > $log_name 2>&1
exit_code="$?"
if [ "$exit_code" -ne "0" ]; then (2>&1 echo "${ERROR} Something got wrong during the execution of $bin_path"); exit -1; fi
echo "${LOG}[$(date)] $bin_name finished"

