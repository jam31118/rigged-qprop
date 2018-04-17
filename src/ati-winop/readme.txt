Example: Hydrogen, lambda=535 nm, I=2e13 W/cm^-2, linear polarization (for comparison to the t-SURFF result)

Quickstart:
make clean
make hydrogen_im
./hydrogen_im
make hydrogen_re
./hydrogen_re
make winop
./winop


Remarks:
"make clean" leads to a fresh compilation of the qprop library. If you run into problems during compilation or execution try this first.

"make hydrogen_im" builds the program for calculationg the initial state by imaginary propagation

"make hydrogen_re" builds the program for the real time propagtion.

"make winop" builds the program for computing the spectrum from the data produced during the real time propagation.

Gnuplot scripts for plotting the results can be found in ../ati-tsurff .
