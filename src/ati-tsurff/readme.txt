Example: Hydrogen, lambda=535 nm, I=2e13 W/cm^-2, linear polarization

Quickstart:
make clean
make hydrogen_im
./hydrogen_im
make hydrogen_re
./hydrogen_re
make eval-tsurff
./eval-tsurff
gnuplot plot-total-spectrum.gp
gnuplot plot-partial-spectra.gp

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

Gnuplot scripts for plotting the total spectrum and the partial spectra for ell=3 and ell=4 are provided.
If the results for the example in ../ati-winop exist these will be plotted as well.
