reset

set term x11

# set term epslatex size 20cm,8cm standalone
# set output "tsurff-mom-res-2000nm-6cyc-1e14-vector.tex"

set term png size 1600,1200
set output "tsurff-mom-res.png"

# set format cb "$%.1t\\times10^{%T}$"
set format cb "$10^{%T}$"

set palette file 'white-blue-red-black-white.gp_palette' using ($1/255):($2/255):($3/255)

set mapping cylindrical
set logscale cb
set pm3d interpolate 4,4 map
# set pm3d map
# ratio of y-axis size to x-axis size
set size ratio 0.33333


set xlabel "parallel momentum (a.u.)"
set ylabel "perpendicular momentum (a.u.)"

# set cbrange [1e-10:1.0]
set xrange [-2:2]
set yrange [:1.5]

# gnuplot expects: theta, z, r
# momentum resolved
splot "< cat -s $(ls -v tsurff-polar*)" every 1:1 u ($3):($4):($2) w pm3d t ""
# energy resolved
# splot "< cat -s $(ls -v tsurff-polar*)" every 1:1 u ($3):($4):($1) w pm3d t ""

set output