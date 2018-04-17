#!/bin/bash

## TO DO
# Move combining program to common places and set environmental variable for qprop or just for this


## Print time
timeOfInitiation=$(date)
echo "[ LOG ] Time of initiation of this script: $timeOfInitiation"


## Parse command-line argument, overiding default configurations
if [ -n "$1" ]; then numOfProcess="$1"; else echo "No first argument \$1, the number of process."; exit 1; fi 

if [ -n "$2" ] && [ "$2" != "." ]; then im_bin="$2"; echo "imaginary propagation binary is set to $im_bin"; fi
if [ -n "$3" ] && [ "$3" != "." ]; then re_bin="$3"; echo "real propagation binary is set to $re_bin"; fi
if [ -n "$4" ] && [ "$4" != "." ]; then eval_bin="$4"; echo "evaluation binary is set to $eval_bin"; fi


## Compile all the binaries
if [[ ! -f "Makefile" && ! -f "makefile" ]]; then echo "No Makefile"; exit 1; fi
echo "[ LOG ] Starting compilation"
make > make.log 2>&1
if [ "$?" -ne 0 ]; then echo "[ERROR] Failed to compile"; exit 1; fi
make eval-tsurff-mpi >> make.log 2>&1
if [ "$?" -ne 0 ]; then echo "[ERROR] Failed to compile"; exit 1; fi
echo "[ LOG ] Completed compilation"


## Basic Configuration
# Simulation program (binary) name list
im_bin=./hydrogen_im
re_bin=./hydrogen_re
eval_bin=./eval-tsurff-mpi

# Parameter Files
initialParamFile="./initial.param"
evalParamFile="./tsurff.param"

# Log files list
im_log="im.log"
re_log="re.log"
eval_log="eval.log"

# MPI setting
numOfProcess=2

# Data file combining program
# [NOTE] (no need at 170926-1642 version)
#combine=./combine_data.sh 

# Data file prefix and suffix
prefix_partial="tsurff-partial"
prefix_polar="tsurff-polar"
suffix=".dat"

# Utility
#tsurffChunkByChunk="../../../common/dispatch.sh"
tsurffChunkByChunk="$QPROP_UTIL/dispatch.sh"
tsurffChunkByChunk_log="./dispatch.log"



## Parse parameter file
expansion_scheme=$(grep -oP 'expansion-scheme long \K(.+)' $evalParamFile)


## Check prerequisites for this code
# Check whether 'mpirun' is present in this system
binName=mpirun
hash $binName 2>/dev/null || { echo >&2 "[ERROR] No $binName found, terminating."; exit 1; }


## Run simulation programs(binaries)
if [ ! -f "$im_bin" ]; then echo "[ERROR] No binary file $im_bin"; exit 1; fi
echo "[ LOG ] Started to run $im_bin"
$im_bin > $im_log 2>&1
if [ "$?" -ne 0 ]; then echo "[ERROR] $im_bin has exited with error code $?, terminating the program"; exit 1; fi
echo "[ LOG ] $im_bin completed with time: $(date)"

if [ ! -f "$re_bin" ]; then echo "[ERROR] No binary file $re_bin"; exit 1; fi
echo "[ LOG ] Started to run $re_bin"
$re_bin > $re_log 2>&1
if [ "$?" -ne 0 ]; then echo "[ERROR] $re_bin has exited with error code $?, terminating the program"; exit 1; fi
echo "[ LOG ] $re_bin completed with time: $(date)"

if [ ! -f "$eval_bin" ]; then echo "[ERROR] No binary file $eval_bin"; exit 1; fi
#echo "[ LOG ] Started to run $eval_bin with $numOfProcess process(es)"
echo "[ LOG ] Started to run $tsurffChunkByChunk with $numOfProcess process(es)"
#mpirun -np $numOfProcess $eval_bin > $eval_log 2>&1
$tsurffChunkByChunk $numOfProcess > $tsurffChunkByChunk_log 2>&1
#if [ "$?" -ne 0 ]; then echo "[ERROR] $eval_bin has exited with error code $?, terminating the program"; exit 1; fi
if [ "$?" -ne 0 ]; then echo "[ERROR] $tsurffChunkByChunk has exited with error code $?, terminating the program"; exit 1; fi
#echo "[ LOG ] $eval_bin completed"
echo "[ LOG ] $tsurffChunkByChunk completed with time: $(date)"


## Combine all data files by multiple processes
# Define function for combining data
function combine_data {
## Command-line argument description:
# $1 -> the number of process used, compulsary
# $2 -> prefix of the data file, compulsary
# $3 -> suffix(file extension specifier) of the data file, default is ".dat"

## Parsing command-line arguments
# Parse the numOfProcess (the number of processes)
if [ -n "$1" ]; then numOfProcess="$1"; 
else
  echo "[ERROR] No argument for \$1 (the number of processes)"
  exit 1
fi

# Parse the prefix of the data filename
if [ -n "$2" ]; then prefix="$2";
else
  echo "[ERROR] No argument for \$2 (the prefix of data filename)"
  exit 1
fi

# Parse the suffix (file extension specifier) of the data filename
if [ -n "$3" ]; then suffix="$3";
else
  echo "[ LOG ] No argument for \$3 (the suffix of data filename)"
  echo "[ LOG ] Falling back to default as \".dat\""
  suffix=".dat"
fi

# Actual concatenation happens here
result=$prefix$suffix
echo -n "" > $result  # make the file empty
for ((id=0; id<$numOfProcess; id++)) {
  obj=$prefix$id$suffix
  if [ -f $obj ]; 
  then 
    #echo "[ LOG ] cat $obj >> $result"; 
    cat $obj >> $result;
  else echo "[ERROR] No file \"$obj\""; exit 1; fi
}

} # function definition ends here


# [NOTE 171002] Replaced by 'dispatch.sh'
# Run function for combining data
#echo "[ LOG ] Started combining data files"
#combine_data $numOfProcess $prefix_polar $suffix
#if [ $expansion_scheme -eq 2 ]; then combine_data $numOfProcess $prefix_partial $suffix; fi
#echo "[ LOG ] Completed combining data files"

echo "[ LOG ] All processes Completed"
timeOfCompletion=$(date)
echo "[ LOG ] Time Of Completion: $timeOfCompletion"

## Implementation using external bash script, invoking another process, which isn't that beautiful.
## Check whether 'bash' is present in this system
#hash bash 2>/dev/null || { echo >&2 "[ERROR] No bash found, terminating."; exit 1; }
## Invoke combining processes
#echo "[ LOG ] Started combining data files"
#bash $combine $numOfProcess $prefix_polar .dat
#if [ "$?" -ne 0 ]; then echo "[ERROR] $combine has exited with error code $?, failed to combine $prefix_polar"; exit 1; fi
#bash $combine $numOfProcess $prefix_partial .dat 
#if [ "$?" -ne 0 ]; then echo "[ERROR] $combine has exited with error code $?, failed to combine $prefix_partial"; exit 1; fi
#echo "[ LOG ] All processes Completed"
