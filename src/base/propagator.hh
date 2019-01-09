#ifndef _PROPAGATOR_HH_
#define _PROPAGATOR_HH_

#include <complex>
#include <cstdlib>

#include "tridiag-common.hh"
//// Define mapping from index_name to integer as indices
//enum index_name { 
//  i_ld, // an index for lower offdiagonal array
//  i_ud, // an index for upper offdiagonal array
//  i_d, // an index for diagonal array
//};


//// Declare functions
void evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis_simple(long l, long sign, long N_rho, double Z, double delta_rho, double delta_t, double *rho_array, double *scalarpot, double *imagpot, std::complex<double> *unitary_tridiags[]);

void evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(long l, long sign, long N_rho, double Z, double delta_rho, double delta_t, double *rho_array, double *scalarpot, double *imagpot, std::complex<double> *diag_unitary, std::complex<double> *lower_offdiag_unitary, std::complex<double> *upper_offdiag_unitary);

void evaluate_diag_V_l(long l, double *scalarpot, double *rho_array, double *imagpot, std::complex<double> *diag_V_l, long N_rho);


#endif // _PROPAGATOR_HH_
