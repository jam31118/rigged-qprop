#ifndef _PROPAGATOR_HH_
#define _PROPAGATOR_HH_

// standard headers
#include <complex>
#include <cstdlib>

// home-made headers
#include "tridiag-common.hh"
#include "matrix.hh"



//// Declare functions

// Generating propagator tridiagonals arrays
void evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis_simple(long l, long sign, long N_rho, double Z, double delta_rho, double delta_t, double *rho_array, double *scalarpot, double *imagpot, std::complex<double> *unitary_tridiags[]);
void evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(long l, long sign, long N_rho, double Z, double delta_rho, double delta_t, double *rho_array, double *scalarpot, double *imagpot, std::complex<double> *diag_unitary, std::complex<double> *lower_offdiag_unitary, std::complex<double> *upper_offdiag_unitary);
void evaluate_diag_V_l(long l, double *scalarpot, double *rho_array, double *imagpot, std::complex<double> *diag_V_l, long N_rho);

// Propagation routine
int crank_nicolson_with_tsurff(
    int index_at_R, double delta_rho, int start_time_index, int num_of_time_steps, 
    std::complex<double> *wf_read, int num_of_wf_lm, 
    std::complex<double> *psi_R_arr, std::complex<double> *dpsi_drho_R_arr, long N_rho, 
    std::complex<double> *tridiags_unitary_stack[], std::complex<double> *tridiags_unitary_inv_stack[], 
    int num_of_steps_to_print_progress, int rank);



#endif // _PROPAGATOR_HH_
