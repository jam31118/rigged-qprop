#ifndef _EVAL_TSURFF_H_
#define _EVAL_TSURFF_H_

#include <memory>
#include <complex>
#include <string>
#include <iostream>

typedef std::complex<double> cplxd;
typedef std::unique_ptr<cplxd[]> cplxd_ptr;

#include <tsurffSpectrum.hh>
#include <potentials.hh>

#include <vecpot.hh>

extern "C" int eval_tsurff(int argc, char **argv);

#endif // _EVAL_TSURFF_H_
