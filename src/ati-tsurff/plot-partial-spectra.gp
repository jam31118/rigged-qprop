reset
set terminal x11

set term png size 1600,1200
set output "partial-spectra.png"

set format y "%.1te%T"

set logscale y

set tmargin at screen 0.95
set bmargin at screen 0.1
set rmargin at screen 0.95
set lmargin at screen 0.1

set yrange [1e-18:1e-2]
set xlabel "electron energy (Hartree)"
set ylabel "differential ionization probability"
set xrange [0:0.72]


plot "< cat -s $(ls -v tsurff-partial*)" u ($1):6 w l t "t-SURFF $\\ell=3$",\
     "< cat -s $(ls -v tsurff-partial*)" u ($1):7 w l t "t-SURFF $\\ell=4$",\
     "../ati-winop/spectrum_0.dat" u 1:($6) w l  t "winop $\\ell=3$",\
     "../ati-winop/spectrum_0.dat" u 1:($7) w l  t "winop $\\ell=4$"

set output
