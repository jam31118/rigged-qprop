#!/bin/bash

# Compilation
make -j4 ppp -C $QPROP_HOME/src/main
make_success="$?"
if [ "$make_success" -ne "0" ]; then echo "[ERROR] during 'make'"; exit -1; fi

if [ -d "ati-ppp" ]; then rm -rf ati-ppp; fi
cp -r ati-before-ppp ati-ppp
cd ati-ppp
$QPROP_HOME/src/main/ppp
ppp_success="$?"
cd ..
if [ "$ppp_success" -ne "0" ]; then exit -1; fi

source ~/py/venv/sci/bin/activate
python ./compare-raw.py
