#include <iostream>
#include <memory>
#include <complex>
#include <exception>
typedef std::complex<double> cplxd;
typedef std::unique_ptr<cplxd[]> cplxd_ptr;

#include <grid.h>
#include <hamop.h>
#include <wavefunction.h>
#include <parameter.hh>
#ifdef HAVE_BOOST
#include <boost/timer.hpp>
#endif
#include <smallHelpers.hh>
#include <powers.hh>
#include <tsurffSpectrum.hh>

#include <vecpot.hh>

// Functions, which determine potentials
#include "potentials.hh"
#include "print-imagpot.hh"

extern "C" int real_prop(int argc, char **argv);

