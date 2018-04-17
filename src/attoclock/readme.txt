Example: Hydrogen, lambda=400 nm, I=1e14 W/cm^-2, circular polarization

Quickstart:
make clean
make hydrogen_im
./hydrogen_im
make hydrogen_re
./hydrogen_re
make eval-tsurff
./eval-tsurff
bash select-omega.sh
gnuplot plot-polar-spectrum.gp


Remarks:
"make clean" leads to a fresh compilation of the qprop library. If you run into problems during compilation or execution try this first.

"make hydrogen_im" builds the program for calculationg the initial state by imaginary propagation

"make hydrogen_re" builds the program for the real time propagtion.

"make eval-tsurff" builds the program for computing the spectrum from the data produced during the real time propagation.

Alternatively a parallel version of the evaluation program can be compiled by "make eval-tsurff-mpi".
It was tested with the openmpi library but using a different mpi library should be possible.
If you want to use a different mpi library or have your openmpi library in a non standard location you need to modify the Makefile.
To start the program with two processes use "mpirun -np 2 eval-tsurff-mpi > log.tsurff".
Every process writes results to its own file so these files need to be joined after the calculation.
For example: cat tsurff-polar0.dat tsurff-polar1.dat  > tsurff-polar.dat
Concatenate the files in the correct order (increasing process number) or you're gonna have a bad time.
Alternativly use "< cat $(ls -v tsurff-polar*.dat)" as the file name in your gnuplot script.

"select-theta.sh" is a bash script that filters only the result for \theta=\pi/2 .

A Gnuplot script for plotting the angle resolved spectrum is provided.
