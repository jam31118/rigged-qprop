#ifndef _EVAL_TSURFF_MPI_H
#define _EVAL_TSURFF_MPI_H

#include <memory>
#include <complex>
#include <string>
#include <iostream>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

typedef std::complex<double> cplxd;
typedef std::unique_ptr<cplxd[]> cplxd_ptr;

#include <tsurffSpectrum.hh>
#include <potentials.hh>

#include <vecpot.hh>

extern "C" int eval_tsurff_mpi(int argc, char **argv);

#endif // _EVAL_TSURFF_MPI_H
