#!/bin/bash
cd ./ati
$QPROP_HOME/bin/imag-prop > im.log 2>&1
$QPROP_HOME/bin/real-prop > re.log 2>&1 &
cd ..

cd ./ati-before-ppp
$QPROP_HOME/bin/imag-prop > im.log 2>&1
$QPROP_HOME/bin/real-prop > re.log 2>&1 &
cd ..
