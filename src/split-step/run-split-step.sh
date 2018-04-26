calc_bin=./hydrogen_re
log_file=./re.log
timesteps_per_iter=10000
total_timesteps=23392
num_of_iter=$(( total_timesteps / timesteps_per_iter + 1 ))
param_file=./propagate.param

printf "num_of_iter: $num_of_iter\n"

echo -n "" > $log_file  # empty the log file

for ((i=0; i<$num_of_iter; i++)) {
  printf "{$i}-th iteration\n";
  start_index=$(( i * timesteps_per_iter ))
  num_index=$timesteps_per_iter
  if [ "$i" -eq "$((num_of_iter - 1))" ]; then num_index=$(( total_timesteps - start_index )); fi
  printf "start_index: $start_index, num_of_timestep: $num_index\n"

  sed -i.bak "s/start-time-index long .*/start-time-index long $start_index/" $param_file
  sed -i.bak "s/num-time-index long .*/num-time-index long $num_index/" $param_file
  grep "^start-time-index" $param_file
  grep "^num-time-index" $param_file

  $calc_bin >> $log_file 2>&1
}

rm $param_file.bak

