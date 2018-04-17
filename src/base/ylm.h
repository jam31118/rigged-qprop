#ifndef ylm_h
#define ylm_h ylm_h
#include<iostream>
#include<cmath>
#include<complex>
#include<factorial.h>

typedef std::complex<double> cplxd;

using namespace std;

cplxd ylm(long l, long m, double theta, double phi);
cplxd ylm2(long l, long m, double theta, double phi);

#endif // ylm_h





