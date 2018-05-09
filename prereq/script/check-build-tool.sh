#!/bin/bash

source $(dirname "$0")/colors.sh

build_tool_commands="gcc g++ make autoconf"
printf "${LOG} Checking for existing build tools: $build_tool_commands\n"
for build_tool_command in $build_tool_commands
do
  path_to_tool_binary=$(command -v $build_tool_command)
  if [ -n "$path_to_tool_binary" ]
  then printf "${OK} $build_tool_command exists at $path_to_tool_binary\n"
  else
    printf "${ERROR} $build_tool_command doesn't exist. Please install it.\n"
    exit -1
  fi
done

