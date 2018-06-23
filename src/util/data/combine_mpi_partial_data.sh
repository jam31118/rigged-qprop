#!/bin/bash


## Source required variables
source $QPROP_HOME/src/util/shell/colors.sh


## Define required functions for combining partial data file produced by each MPI process
function check_and_combine_mpi_data {
  ## Assign function arguments to variables, which are local in function
  local _obj_prefix="$1"
  local _obj_suffix="$2"

  local _obj="$_obj_prefix.$_obj_suffix"

  ## Check whether the combined data is already there
  if [ -f "$_obj" ]; then echo "${LOG} There is $_obj alreadly."
  echo "${LOG} deleting $_obj . . ."; rm "$_obj"
  else echo "${LOG} The file doesn't exists: $_obj."; fi
  
  ## Check whether partial data files exist
  local _partial_data_filename_format="$_obj_prefix"*."$_obj_suffix"
  local _num_of_partial_data_file=$(ls -1 $_partial_data_filename_format 2> /dev/null | wc -l)
  if [ "$_num_of_partial_data_file" -eq "0" ]
  then echo "${LOG} No data file is found: $_partial_data_filename_format will be omitted."; return 0; fi

  ## Combine partial data files
  for f in $_partial_data_filename_format
  do
  	echo "${LOG} Combining \"$f\" to \"$_obj\""
  	cat "$f" >> "$_obj"
  done

  echo "${OK} Combining partial data succeeded for $_partial_data_filename_format"

  return 0
}

