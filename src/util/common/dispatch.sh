#!/bin/bash

## NOTE ##
# This script only works for tsurff calculation
# WINOP version will be implemented later


evalParamFile="./tsurff.param"
evalParamFileCopy="./tsurffpsi.param.copy"


## Parse parameter file
expansion_scheme=$(grep -oP 'expansion-scheme long \K(.+)' $evalParamFile)

## Determine file type list to concathenate
concatList="tsurff-polar"
if [ $expansion_scheme -eq 2 ]; then concatList=$concatList" tsurff-partial"; fi


## Parse command-line arguments
numOfProcess=$1
if [ -z "$numOfProcess" ]; then echo "[ERROR] No numOfProcess input"; exit 1; fi

eval_bin="./eval-tsurff-mpi"
if [ -n "$2" ]; then eval_bin="$2"; fi

#requiredTotalMemCalculator="../../../common/getTotalSizeOf_tsurff.sh"
requiredTotalMemCalculator="$QPROP_UTIL/getTotalSizeOf_tsurff.sh"
if [ -n "$3" ]; then requiredTotalMemCalculator="$3"; fi


## Make resulting data file empty
suffix=".dat"
for prefix in $concatList
do
  #result_data_polar="./tsurff-polar.dat"
  #echo -n "" > $result_data_polar

  result_data="./$prefix$suffix"
  echo -n "" > $result_data
done 


# Get parameter from parameter file
num_theta_surff=$(grep -oP 'num-theta-surff long \K(.+)' $evalParamFile)

# Memory threshold that a program (including all processes of MPI) shouldn't use beyond
possibleMem=3500000000 

# Total amount of memory the program will use (estimated minimum, possibly more)
requiredTotalMem=$(bash $requiredTotalMemCalculator)

## Required number of calculations = ceil( Required Total Memory / Possible(i.e. available) Memory )
## .. where 'ceil' mean 'round up', e.g. ceil(3.2) == 3
## .. and 'ceil(A/B)' is equivalent to '(A+B-1)/A' where '/' denotes integer division
#requiredNumOfCalcMinimum=$(echo "($requiredTotalMem + $possibleMem - 1) / $possibleMem" | bc -l)
requiredNumOfCalcMinimum=$(( ($requiredTotalMem + $possibleMem - 1) / $possibleMem ))

## Get the number of theta-values which should be calculated in each calculation
## .. number of Theta Per Calculation = floor( total number of Theta / required Number of Calculation )
#numThetaPerCalc=$(( ($num_theta_surff + $requiredNumOfCalcMinimum - 1) / $requiredNumOfCalcMinimum  ))
numThetaPerCalc=$(( $num_theta_surff / $requiredNumOfCalcMinimum  ))
if [ "$numThetaPerCalc" -eq 0 ]; then (echo "[ERROR] Cannot divide workload anymore w.r.t theta"); exit 1; fi

requiredNumOfCalc=$(( ($num_theta_surff + $numThetaPerCalc - 1) / $numThetaPerCalc ))


echo numThetaPerCalc = $numThetaPerCalc
echo num_theta_surff = $num_theta_surff
echo requiredTotalMem = $requiredTotalMem
echo requiredNumOfCalcMinimum = $requiredNumOfCalcMinimum
echo requiredNumOfCalc = $requiredNumOfCalc 

temp_data_dir="temp_data"
if [ ! -d "$temp_data_dir" ]; then mkdir $temp_data_dir; fi

min_theta_index=0
num_theta_partial=0
curIndexOfCalc=0
while [ "$curIndexOfCalc" -lt "$requiredNumOfCalc" ]
do
  #echo curIndexOfCalc = $curIndexOfCalc
  min_theta_index=$(echo "$curIndexOfCalc * $numThetaPerCalc" | bc -l)
  let num_theta_partial=$numThetaPerCalc
  
  let next_min_theta_index=($min_theta_index + $numThetaPerCalc)
  if [ "$next_min_theta_index" -gt "$num_theta_surff" ];
  then
    num_theta_partial=$(echo "$num_theta_surff - $min_theta_index" | bc -l)
  fi
  echo min_theta_index = $min_theta_index
  echo num_theta_partial = $num_theta_partial

#  mv $evalParamFile $evalParamFileCopy

  # Backup original parameter file
  cp $evalParamFile $evalParamFileCopy  

  # Modify parameter file(s)
  cat $evalParamFile \
    | sed "s/min-theta-index long .*/min-theta-index long $min_theta_index/g"  \
    | sed "s/num-theta-partial long .*/num-theta-partial long $num_theta_partial/g" \
    > $evalParamFile$curIndexOfCalc

  cp $evalParamFile$curIndexOfCalc $evalParamFile
  
  echo "[ LOG ] $eval_bin started with time: $(date)"
  mpirun -np $numOfProcess $eval_bin > eval.log.$curIndexOfCalc 2>&1
  if [ "$?" -ne 0 ]; then echo "[ERROR] $eval_bin has exited with error code $?, terminating the program"; exit 1; fi
  echo "[ LOG ] $eval_bin completed with time: $(date)"
  
  temp_data_subdir="$temp_data_dir/$curIndexOfCalc"
  if [ ! -d "$temp_data_subdir" ]; then mkdir $temp_data_subdir; fi
  #mv ./tsurff-polar*.dat ./tsurff-partial*.dat "$temp_data_subdir" # no tsurff-partial file in expansion-scheme = 1
  mv ./tsurff-polar*.dat "$temp_data_subdir"
  if [ $expansion_scheme -eq 2 ]; then mv ./tsurff-partial*.dat "$temp_data_subdir"; fi

  # Only for expansion-scheme = 1 which doesn't generate tsurff-partial
  # [171017 NOTE] -> fixed to include also for expansion-scheme == 2

  # Actual concatenation happens here
  suffix=".dat"
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
  
  #cat "$subresult" >> "$result_data_polar"
  #echo "[ LOG ] cat $subresult >> $result_data_polar"

  curIndexOfCalc=$(($curIndexOfCalc + 1))
done




