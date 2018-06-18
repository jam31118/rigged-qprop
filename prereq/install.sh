#!/bin/bash

## Source required scripts
source $(dirname "$0")/script/colors.sh


## Check environment variables setting
if [ ! -d "$QPROP_HOME" ]
then
  (>&2 echo -e "${ERROR} Please set \$QPROP_HOME to valid path where source code resides")
  (>&2 echo -e "${LOG} Now, \$QPROP_HOME is set as: $QPROP_HOME")
  (>&2 echo -e "${LOG} If you have set a variable but it doesn't seem to appear in this script, ")
  (>&2 echo -e "${LOG} check if the variable has been exported using 'export' command.")
  exit -1
fi

if [ ! -d "$QPROP_DEP_DIR" ]
then
  (>&2 echo -e "${ERROR} Please set \$QPROP_DEP_DIR to valid path where you want to install dependencies of QPROP")
  exit -1
fi


## Check whether required source code are available
url_gsl="http://ftp-stud.hs-esslingen.de/pub/Mirrors/ftp.gnu.org/gsl/gsl-2.4.tar.gz"
url_openmpi="https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.0.tar.gz"
url_boost="https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz"
for url in $url_gsl $url_openmpi $url_boost
do
  wget --spider $url
  wget_exit_code="$?"
  if [ "$wget_exit_code" -ne "0" ]
  then
    (>&2 echo "${ERROR} Cannot find requested webpage: $url")
    (>&2 echo "${ERROR} Please check if the source code in the webpage is available in the webpage")
  else echo "${LOG} Found source code at $url"
  fi
done

## Prepare relavant directories
BASE_DIR="$QPROP_DEP_DIR"

SCRIPT_DIR="$QPROP_HOME/prereq/script"
if [ ! -d "$QPROP_HOME" ]; then (>&2 echo -e "${ERROR} Scripts doens't exist"); exit -1; fi

LOG_DIR="$BASE_DIR/log"
if [ ! -d "$LOG_DIR" ]; then mkdir -p $LOG_DIR; fi


## Logging relavant directories
echo -e "${LOG} BASE_DIR: $BASE_DIR"
echo -e "${LOG} SCRIPT_DIR: $SCRIPT_DIR"
echo -e "${LOG} LOG_DIR: $LOG_DIR"


## Check build tools
bash $SCRIPT_DIR/check-build-tool.sh
if [ "$?" -ne "0" ]; then (>&2 printf "${ERROR} Failed to find required build tools.\n"); exit -1; fi


## Install
for program in gsl openmpi boost
do
  script_path="$SCRIPT_DIR/install-$program.sh"
  if [ ! -f "$script_path" ]; then (>&2 echo -e "${ERROR} Script \'$script_path\' doens't exist"); exit -1; fi
  log_file_path="$LOG_DIR/$program.log"
  echo -e "${LOG} Installing $program . . . "
  echo -e "${LOG} For detail information, refer to $log_file_path"
  bash $script_path $BASE_DIR >&1 2>&1 | tee $log_file_path
  if [ "${PIPESTATUS[0]}" -ne "0" ]
  then (>&2 echo -e "${ERROR} Failed to install $program"); exit -1; 
  else printf "${OK} Installing $program succeeded.\n"
  fi
done

rm -rf $LOG_DIR

exit 0
