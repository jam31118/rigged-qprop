
#!/bin/bash

## Check environment variables setting
if [ ! -d "$QPROP_HOME" ]
then
  (>&2 echo "[ERROR] Please set \$QPROP_HOME to valid path where source code resides")
  exit -1
fi

if [ ! -d "$QPROP_DEP_DIR" ]
then
  (>&2 echo "[ERROR] Please set \$QPROP_DEP_DIR to valid path where you want to install dependencies of QPROP")
  exit -1
fi


BASE_DIR="$QPROP_DEP_DIR"

SCRIPT_DIR="$QPROP_HOME/prereq/script"
if [ ! -d "$QPROP_HOME" ]; then (>&2 echo "[ERROR] Scripts doens't exist"); exit -1; fi

LOG_DIR="$BASE_DIR/log"
if [ ! -d "$LOG_DIR" ]; then mkdir -p $LOG_DIR; fi


echo "[ LOG ] BASE_DIR: $BASE_DIR"
echo "[ LOG ] SCRIPT_DIR: $SCRIPT_DIR"
echo "[ LOG ] LOG_DIR: $LOG_DIR"


for program in gsl openmpi boost
do
  script_path="$SCRIPT_DIR/install-$program.sh"
  if [ ! -f "$script_path" ]; then (>&2 echo "[ERROR] Script \'$script_path\' doens't exist"); exit -1; fi
  log_file_path="$LOG_DIR/$program.log"
  echo "[ LOG ] Installing $program . . . "
  echo "[ LOG ] For detail information, refer to $log_file_path"
  bash $script_path $BASE_DIR >&1 2>&1 | tee $log_file_path
  if [ "$?" -ne "0" ]; then (>&2 echo "[ERROR] Failed to install $program"); exit -1; fi
done

rm -rf $LOG_DIR

exit 0
