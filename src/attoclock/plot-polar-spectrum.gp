reset
set term x11

set term png size 1600,1200
set output "polar-spectrum.png"

set palette file 'white-blue-red-black-white.gp_palette' using ($1/255):($2/255):($3/255)

# Heaviside function
theta(x)=(x<0)?0.0:1.0

set xlabel "momentum $x$-direction"
set ylabel "momentum $y$-direction"

set mapping cylindrical
set pm3d map

# set logscale cb

set xrange [-1.0:1.0]
set yrange [-1.0:1.0]
set tmargin  at screen 0.95
set bmargin  at screen 0.1
set lmargin at screen 0.1
set rmargin at screen 0.9

# gnuplot expects: theta, z, r
splot "tsurff-polar.dat" u 4:($5)*(theta(1.0-$2)):2 w pm3d t ""

set output