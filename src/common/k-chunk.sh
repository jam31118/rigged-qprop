#!/bin/bash

source $QPROP_UTIL/param.sh

### Parse command-line arguments
if [ -n "$1" ]; then numOfProcess="$1"
else (>&2 echo "[ERROR] No numOfProcess input"); exit 1; fi

if [ -n "$2" ]; then memLimit_GB="$2"
else (>&2 echo "[ERROR] No memory limit given"); exit 1; fi
#defaultMemoryLimit=00


### Check required utilities and files
evalParamFile="./tsurff.param"
requiredTotalMemCalculator="$QPROP_UTIL/getTotalSizeOf_tsurff.sh"
eval_bin="./eval-tsurff-mpi"
mpi_bin="mpirun"

util_list="$evalParamFile $requiredTotalMemCalculator"
for util in $util_list
do 
  if [ ! -e $util ]; then (>&2 echo "[ERROR] Following utility doesn't exist: $util"); exit 1; fi
done


## Other configurations
evalParamFileCopy="./tsurffpsi.param.copy"
temp_data_dir="temp_data"
if [ ! -d "$temp_data_dir" ]
  then mkdir $temp_data_dir
  else rm -rf $temp_data_dir
fi


## Parse parameter file
expansion_scheme=$(get_param expansion-scheme $evalParamFile)
if [ $? -ne 0 ]; then (>&2 echo "[ERROR] Failed to parse parameter"); exit 1; fi
num_k_surff=$(get_param num-k-surff $evalParamFile)
if [ $? -ne 0 ]; then (>&2 echo "[ERROR] Failed to parse parameter"); exit 1; fi


## Determine file type list to concathenate
concatList="tsurff-polar"
if [ $expansion_scheme -eq 2 ]; then concatList=$concatList" tsurff-partial"; fi


## Make resulting data file empty
suffix=".dat"
for prefix in $concatList
do
  result_data="./$prefix$suffix"
  echo -n "" > $result_data
done 



# Memory threshold that a program (including all processes of MPI) shouldn't use beyond
possibleMem_float=$(echo "$memLimit_GB * 1024 * 1024 * 1024" | bc -l)
possibleMem=$(LC_NUMERIC="en_US.UTF-8" printf "%.0f" $possibleMem_float)
#possibleMem=$(( $memLimit_GB * 1024 * 1024 )) # in unit of Bytes

# Total amount of memory the program will use (estimated minimum, possibly more)
requiredTotalMem_min=$(bash $requiredTotalMemCalculator)
buffer_factor="1.5"
requiredTotalMem_float=$(echo "$requiredTotalMem_min * $buffer_factor" | bc -l)
requiredTotalMem=$(LC_NUMERIC="en_US.UTF-8" printf "%.0f" $requiredTotalMem_float)


## Required number of calculations = ceil( Required Total Memory / Possible(i.e. available) Memory )
## .. where 'ceil' mean 'round up', e.g. ceil(3.2) == 3
## .. and 'ceil(A/B)' is equivalent to '(A+B-1)/A' where '/' denotes integer division
#requiredNumOfCalcMinimum=$(echo "($requiredTotalMem + $possibleMem - 1) / $possibleMem" | bc -l)
requiredNumOfCalcMinimum=$(( ($requiredTotalMem + $possibleMem - 1) / $possibleMem ))


## Get the number of k-values which should be calculated in each calculation
## .. number of k-values Per Calculation = floor( total number of k-values / required Number of Calculation )
num_k_per_calc=$(( $num_k_surff / $requiredNumOfCalcMinimum  ))
if [ "$num_k_per_calc" -eq 0 ]; then (>&2 echo "[ERROR] Cannot divide workload anymore w.r.t k"); exit 1; fi

requiredNumOfCalc=$(( ($num_k_surff + $num_k_per_calc - 1) / $num_k_per_calc ))


## Logging
echo num_k_per_calc = $num_k_per_calc
echo num_k_surff = $num_k_surff
memLimit_MB=$(echo "$memLimit_GB * 1024" | bc -l)
memLimit_MB_int=$(LC_NUMERIC="en_US.UTF-8" printf "%.0f" $memLimit_MB)
if [ $memLimit_MB_int -eq 0 ]; then memLimit_MB_int=$memLimit_MB; fi
echo "memLimit_MB = $memLimit_MB_int MB"
requiredTotalMem_MB=$(echo "$requiredTotalMem / 1024 / 1024" | bc -l)
requiredTotalMem_MB_int=$(LC_NUMERIC="en_US.UTF-8" printf "%.0f" $requiredTotalMem_MB)
if [ $requiredTotalMem_MB_int -eq 0 ]; then requiredTotalMem_MB_int=$requiredTotalMem_MB; fi
echo "requiredTotalMem = $requiredTotalMem_MB_int MB"
echo requiredNumOfCalcMinimum = $requiredNumOfCalcMinimum
echo requiredNumOfCalc = $requiredNumOfCalc 


## Backup original parameter file
cp $evalParamFile $evalParamFileCopy  


## Check and add for temporary parameters: min-k-index, num-k-partial
temp_param_list="min-k-index num-k-partial"
for temp_param in $temp_param_list
do
  num_of_temp_param=$(grep $temp_param $evalParamFile | wc -l)
  if [ $num_of_temp_param -eq 0 ]; then echo "$temp_param long 0" >> $evalParamFile; fi
done


## Start looping on required repeatition of calculation
min_k_index=0
num_k_partial=0
curIndexOfCalc=0
while [ "$curIndexOfCalc" -lt "$requiredNumOfCalc" ]
do
  min_k_index=$(echo "$curIndexOfCalc * $num_k_per_calc" | bc -l)
  let num_k_partial=$num_k_per_calc
  
  let next_min_k_index=($min_k_index + $num_k_per_calc)
  if [ "$next_min_k_index" -gt "$num_k_surff" ];
  then
    num_k_partial=$(echo "$num_k_surff - $min_k_index" | bc -l)
  fi
  echo min_k_index = $min_k_index
  echo num_k_partial = $num_k_partial


  ## Modify parameter file(s)
  temp_param_file="$evalParamFile.$curIndexOfCalc"

  cat $evalParamFile \
    | sed "s/min-k-index long .*/min-k-index long $min_k_index/g"  \
    | sed "s/num-k-partial long .*/num-k-partial long $num_k_partial/g" \
    > $temp_param_file

  cp $temp_param_file $evalParamFile


  ## Run calculations
  echo "[ LOG ] $eval_bin started with time: $(date)"
  temp_eval_log_file="eval.log.$curIndexOfCalc"
  $mpi_bin -np $numOfProcess $eval_bin > $temp_eval_log_file 2>&1
  status_code=$?
  if [ "$status_code" -ne 0 ]; then echo "[ERROR] $eval_bin has exited with error code $status_code, terminating the program"; exit 1; fi
  echo "[ LOG ] $eval_bin completed with time: $(date)"

  temp_data_subdir="$temp_data_dir/$curIndexOfCalc"
  if [ ! -d "$temp_data_subdir" ]; then mkdir -p $temp_data_subdir; fi
  
  suffix=".dat"
  mv $temp_param_file $temp_eval_log_file "$temp_data_subdir"
  for prefix in $concatList
  do
	  for ((id=0; id<$numOfProcess; id++)) { 
      mv ./$prefix$id$suffix "$temp_data_subdir" 
    }
  done

  # Actual concatenation happens here
  for prefix in $concatList
  do
	  subresult="$temp_data_subdir/$prefix$suffix"
	  echo -n "" > $subresult  # make the file empty
	  for ((id=0; id<$numOfProcess; id++)) {
	    obj="$temp_data_subdir/$prefix$id$suffix"
	    if [ -f $obj ]; 
	    then 
	      echo "[ LOG ] cat $obj >> $subresult" 
	      cat $obj >> $subresult
	    else echo "[ERROR] No file \"$obj\""; exit 1; fi
	  }

	  result_data="./$prefix$suffix"

	  cat "$subresult" >> "$result_data"
	  echo "[ LOG ] cat $subresult >> $result_data"
  done

  ## Remove temporary files
  #rm $temp_param_file
  
  ## Required for looping to end
  curIndexOfCalc=$(($curIndexOfCalc + 1))
done


## Restore original parameter file
mv $evalParamFileCopy $evalParamFile
