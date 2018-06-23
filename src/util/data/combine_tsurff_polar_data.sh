#!/bin/bash

source $QPROP_HOME/src/util/data/combine_mpi_partial_data.sh

obj_prefix="tsurff-polar"
obj_suffix="dat"

check_and_combine_mpi_data $obj_prefix $obj_suffix
filename_format="$obj_prefix*.$obj_suffix"
if [ "$?" -ne "0" ]; then (>&2 echo "${ERROR} Failed to combine partial data for $filename_format"); fi
#else echo "${OK} Combining partial data succeeded for $filename_format"; fi

