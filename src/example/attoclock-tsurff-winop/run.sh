#!/bin/bash


## Step 1
$QPROP_HOME/bin/imag-prop
$QPROP_HOME/bin/real-prop


## Step 2 (select either one or both)

# for tsurff
$QPROP_HOME/bin/ppp
$QPROP_HOME/bin/eval-tsurff
cp tsurff-polar0.dat tsurff-polar.dat
cp tsurff-partial0.dat tsurff-partial.dat

# for winop
$QPROP_HOME/bin/winop
cp spectrum_0.dat spectrum.dat
cp spectrum_polar0.dat spectrum_polar.dat

