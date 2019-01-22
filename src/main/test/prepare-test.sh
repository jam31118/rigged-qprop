#!/bin/bash

source $QPROP_HOME/src/util/shell/colors.sh

make -j4 -C $QPROP_HOME/src/main
exit_code_make=$?
if [ "$exit_code_make" -eq "0" ]
then echo "${LOG} Compliation done.";
else (2>&1 echo "${ERROR} Compilation failed\n"); exit -1;
fi

cd ./ati
$QPROP_HOME/bin/imag-prop > im.log 2>&1
$QPROP_HOME/bin/real-prop > re.log 2>&1 &
cd ..

cd ./ati-before-ppp
$QPROP_HOME/bin/imag-prop > im.log 2>&1
$QPROP_HOME/bin/real-prop > re.log 2>&1 &
cd ..
