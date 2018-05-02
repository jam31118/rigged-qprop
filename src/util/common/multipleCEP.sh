#!/bin/bash
echo "Started $0 with time: $(date)"


## Parse Command-line arguments
if [ -n "$1" ]; then numOfProcess="$1"; else echo "No \$1, the number of processes of mpi"; exit 1; fi
if [ -n "$2" ]; then numOfCEP="$2"; else echo "No \$2, the number of CEP values"; exit 1; fi


## Configurations
propParamFile="./propagate.param"
propParamFileCopy="./propagate.param.copy"


pi=3.141592653589793


if [ "$numOfCEP" -lt 2 ]; then echo "[ERROR] Number of CEP should be two or more"; exit 1; fi
#upperLimitOfCEP=$(echo "2*$pi" | bc -l)
upperLimitOfCEP=$(echo "$pi" | bc -l)
CEP_interval=$(echo "$upperLimitOfCEP / ($numOfCEP - 1)" | bc -l)


indexOfCEP=0
while [ "$indexOfCEP" -lt "$numOfCEP" ]
do
  ## Determine CEP value
  CEP=$(echo "$indexOfCEP * $CEP_interval" | bc -l)
  echo "Current CEP: $CEP"


  ## Requried directory
  data_subdir="cep_data/$indexOfCEP"
  if [ ! -d "$data_subdir" ]; then mkdir -p $data_subdir; fi


  ## Construct propagate parameter file with backup
  cp $propParamFile $propParamFileCopy 

  cat $propParamFile \
    | sed "s/cep double .*/cep double $CEP/g" \
    > "$propParamFile$indexOfCEP"

  cp "$propParamFile$indexOfCEP" $propParamFile 
  cp $propParamFile "$data_subdir/$propParamFile"
  mv "$propParamFile$indexOfCEP" "$data_subdir/$propParamFile$indexOfCEP"

  ## Run calculation for specified CEP 
  printf "indexOfCEP: $indexOfCEP/$numOfCEP\nStarted with time: $(date)\n"
  bash ../../../common/runall.sh $numOfProcess > runall.log 2>&1 
  if [ "$?" -ne 0 ]; then echo "[ERROR] runall.sh exited abnormally"; exit 1; fi
  printf "indexOfCEP: $indexOfCEP/$numOfCEP\nEnded with time: $(date)\n"
  printf "\n"

  ## Move resulting data 
  cp *.raw *.dat *.log *.log.* *.param $data_subdir


  ## For looping
  indexOfCEP=$(($indexOfCEP + 1))
done


echo "All processes of $0 ended with time: $(date)"


exit 0;
