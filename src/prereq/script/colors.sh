#!/bin/bash

#BD="\e[1m"
#R="\033[1;31m"
#G="\033[1;32m"
#NC="\033[0m"

supported_num_of_colors=$(tput colors)
required_num_of_colors=16

if test -t 1 && test -n "$supported_num_of_colors" && test $supported_num_of_colors -ge $required_num_of_colors
then

  BD="$(tput bold)"
  R="$(tput setaf 9)"
  G="$(tput setaf 10)"
  Y="$(tput setaf 11)"
  B="$(tput setaf 12)"
  NC="$(tput sgr0)"

fi

OK="${G}${BD}[ O K ]${NC}"
LOG="${B}${BD}[ LOG ]${NC}"
ERROR="${R}${BD}[ERROR]${NC}"

#echo "$OK"
#echo "$LOG"
#echo "$ERROR"

