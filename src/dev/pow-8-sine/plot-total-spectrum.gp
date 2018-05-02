reset
set term x11

set term png size 1600,1200
set output "total-spectrum.png"

set logscale y
set format y "%.te%T"

set xlabel "electron energy (Hartree)"
set ylabel "differential ionization probability"

set tmargin at screen 0.95
set bmargin at screen 0.1
set rmargin at screen 0.95
set lmargin at screen 0.1

set xrange [0:0.72]

plot\
     "< cat -s $(ls -v tsurff-partial*)" u ($1):($18) w l  t "t-SURFF"

set output
