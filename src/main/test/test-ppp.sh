#!/bin/bash

# Compilation
make ppp -C $QPROP_HOME/src/main

if [ -d "ati-ppp" ]; then rm -rf ati-ppp; fi
cp -r ati-before-ppp ati-ppp
cd ati-ppp
$QPROP_HOME/src/main/ppp
cd ..

source ~/py/venv/sci/bin/activate
python ./compare-raw.py
