#!/bin/bash
# get rid of old result
rm tsurff-polar.dat
# print only the lines for theta=pi/2 and blank lines between the data blocks
awk '$3=="1.5707963267948966" || $0=="" {print $0}' $(ls -v tsurff-polar*) > tsurff-polar.dat
# erase superfluous blank lines
cat -s tsurff-polar.dat > temp
# copy line for phi=0 to the end of a data block
awk '$4==0 { line=$0 }; $0=="" { $0 = $0 line "\n" }; { print $0 }' temp > tsurff-polar.dat
rm temp
