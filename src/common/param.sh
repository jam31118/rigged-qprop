#!/bin/bash

LIGHTRED="\033[1;31m"
NOCOLOR="\033[0m"

## Parameter file parser, specific to Qprop 2.0
function get_param {
  ## Check arguments
  if [ -n "$1" ]; then (>&2 echo "[ LOG ] Given parameter name: $1"); parameter_name="$1"
  else (>&2 echo "[ERROR] No parameter name given"); exit 1; fi
  if [ -n "$2" ]; then (>&2 echo "[ LOG ] Given parameter file name: $2"); param_file_name="$2"
  else (>&2 echo "[ LOG ] No parameter file name given, searching all files with .param extension"); param_file_name="./*.param"; fi
  #else (>&2 echo "[ERROR] No parameter file name given"); exit 1; fi
  
  let numOfMatchedParamEntry=0
  for param_file in $param_file_name
  do
    parsing_result="$(grep -oP "$parameter_name\s+[a-zA-Z0-9-]+\s+\K(.+)" $param_file)"
    #(>&2 echo $parsing_result)
    if [ -n "$parsing_result" ]; then parsed_param="$parsing_result"; let numOfMatchedParamEntry=$numOfMatchedParamEntry+1; fi
    #(>&2 echo $numOfMatchedParamEntry)
  done

  if [ "$numOfMatchedParamEntry" -lt "0" ]; then (>&2 echo -e "${LIGHTRED}[ERROR]${NOCOLOR} Unexpected exception: the number of matching parameter entry is zero"); exit 1; fi
  if [ "$numOfMatchedParamEntry" -eq "0" ]; then (>&2 echo -e "${LIGHTRED}[ERROR]${NOCOLOR} No matching parameter entry in $param_file_name with parameter name: $parameter_name"); exit 1; fi
  if [ "$numOfMatchedParamEntry" -gt "1" ]; then (>&2 echo -e "${LIGHTRED}[ERROR]${NOCOLOR} Multiple entries of same parameter name: $parameter_name"); exit 1; fi
  if [ "$numOfMatchedParamEntry" -eq "1" ]; then (>&2 echo "[ LOG ] Successfully parsed parameter $parameter_name"); fi
  #if [ "$?" -ne 0 ]; then (>&2 echo "[ERROR] Failed to parse $parameter_name from $param_file_name"); exit 1; fi
  echo $parsed_param
}


#result=$(get_param expansion-schem tsurff.param)
#if [ $? -ne 0 ]; then (>&2 echo -e "${LIGHTRED}[ERROR]${NOCOLOR} Failed to parse parameter file")
#else echo result: $result; fi


