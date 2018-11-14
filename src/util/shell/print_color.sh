#!/bin/bash

min_index=0
max_index=16

if [ -n "$1" ]; then max_index="$1"; fi
if [ -n "$2" ]; then min_index="$2"; fi

normal_color=$(tput sgr0)

for ((i=$min_index; i<$max_index; i++))
do 
  color=$(tput setaf $i)
  echo "${color}this color is given by \$(tput setaf $i)${normal_color}"
done

