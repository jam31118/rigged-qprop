#!/bin/bash

binary_file_name="real-prop"
log_file_name="re.log"
column_index_of_dirpath_from_pwdx=2

num_of_line_to_show=3

for p in $(pgrep $binary_file_name)
do
  dirpath=$(pwdx $p | cut -d ' ' -f $column_index_of_dirpath_from_pwdx)
  echo "The log file '$log_file_name' of '$binary_file_name' (pid: $p) at '$dirpath'"
  tail $dirpath/$log_file_name -n $num_of_line_to_show
  echo ""
done

