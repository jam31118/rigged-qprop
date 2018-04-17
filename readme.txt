This is a Qprop release for calculating electron spectra of single electron systems in intense laser fields using the t-SURFF approximation.

Several examples are provided in the src/ sub directory.

Prerequisites for using the package:
-GNU Scientific Library (gsl versions 1.16 and 2.0 have been tested)
-GNU C++ compiler supporting the C++11 standard (Versions 4.7.2, 4.8.4 and 5.3.0 have been tested)

Optional:
-Boost Timer (measure execution times; enable with compiler flag -DHAVE_BOOST)
-Open MPI (parallelized version of spectrum evaluation)


If these libraries are in non standard locations the Makefiles in the example directories and in src/base have to be modified.
