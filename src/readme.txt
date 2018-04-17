The subdirectory base/ contains the source files of the qprop library and a Makefile to build the library.
If you plan to use a different compiler than g++ you need to modify this Makefile.

The example subdirectories:

1. ati-tsurff/
This is the least time consuming example where you calculate a spectrum for above threshold ionization of Hydrogen (lambda=535 nm, I=2e13 W/cm^-2, linear polarization).

2. ati-winop/
Compare the spectrum of 1. to the window operator result.

3. large-clubs/
Calculate  a much more demanding  angle resolved spectrum for Hydrogen.
In this case a large number of photons is involved to create an impressive club structure in the electron spectrum (lambda=2000 nm, I=1e14 W/cm^-2, linear polarization)

4. attoclock/
Here an angle resolved electron spectrum for the ionization of Hydrogen in a circularly polarized two cycle laser pulse is calculated (lambda=400 nm, I=1e14 W/cm^-2, circular polarization).

5. pow-8-sine/
Like ati-tsurff but for a sin^8 pulse envelope. Besides total and partial energy PES a momentum-resolved PES is generated with plot-polar-spectrum.gp.
