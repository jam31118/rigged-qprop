#include "propagator.hh"

void evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(long l, long sign, long N_rho, double Z, double delta_rho, double delta_t, double *rho_array, double *scalarpot, double *imagpot, std::complex<double> *diag_unitary, std::complex<double> *lower_offdiag_unitary, std::complex<double> *upper_offdiag_unitary) {
  /* # Arguments
   * `scalarpot`, `imagpot`, `diag_unitary`: each of length `N_rho`
   * `lower_offdiag_unitary`, `upper_offdiag_unitary`: each of length `N_rho-1`
   * */

  //// Declare required variables
//  std::complex<double> one_i = std::complex<double>(0,1.0);
  std::complex<double> coef = std::complex<double>(0,sign * delta_t * 0.5);
  double h = delta_rho;
  std::complex<double> c1, c2;

  //// Evaulate `diag_V_l`
  std::complex<double> *diag_V_l;
  diag_V_l = new std::complex<double>[N_rho];
  evaluate_diag_V_l(l, scalarpot, rho_array, imagpot, diag_V_l, N_rho);

  std::complex<double> *p_diag_V_l;

  //// Evaulate `diag_unitary`
  c1 = -5.0/3.0 + coef * (-2.0/(h*h));
  std::complex<double> *p_diag_unitary = diag_unitary,
    *p_diag_unitary_max = diag_unitary + N_rho;
  p_diag_V_l = diag_V_l;
  for (; p_diag_unitary < p_diag_unitary_max; ++p_diag_unitary, ++p_diag_V_l) {
    *p_diag_unitary = c1 - 5.0/3.0 * coef * (*p_diag_V_l);
  }
  if (l==0) { diag_unitary[0] += Z*h/(12.0-10.0*Z*h) * (-1.0/3.0 * (1.0+coef*diag_V_l[0]) + 2.0/(h*h)*coef); }

  //// Evaluate `lower_` and `upper_` `offdiag_unitary`
  std::complex<double> *w = new std::complex<double>[N_rho];
  c2 = -1.0/6.0 + coef / (h*h);
  std::complex<double> coef_over_6 = coef / 6.0;
  std::complex<double> *p_w = w, *p_w_max = w + N_rho;
  p_diag_V_l = diag_V_l;
  for (; p_w < p_w_max; ++p_w, ++p_diag_V_l) {
    *p_w = c2 - coef_over_6 * *p_diag_V_l;
  }
  long i;
  lower_offdiag_unitary[0] = w[0];
  for (i=1; i<N_rho-1; i++) {
    lower_offdiag_unitary[i] = w[i];
    upper_offdiag_unitary[i-1] = w[i];
  }
  upper_offdiag_unitary[N_rho-2] = w[N_rho-1];
  // [NOTE] note that `lower_offdiag_unitary` and `upper_offdiag_unitary` is of length `N_rho-`
  
  //// Free temporary array
  delete [] w;
  delete [] diag_V_l;
}

void evaluate_diag_V_l(long l, double *scalarpot, double *rho_array, double *imagpot, std::complex<double> *diag_V_l, long N_rho) {
  // Prepare array pointers for evaulating `diag_V_l`
  double *p_scalarpot = scalarpot,
         *p_rho_array = rho_array,
         *p_imagpot = imagpot;
  std::complex<double> *p_diag_V_l = diag_V_l, *p_diag_V_l_max = diag_V_l + N_rho;
  
  // Prepare temporary variable
  double rho;

  // Evaulate constant first
  double coef = 0.5 * l * (l+1);

  // Evaulate `diag_V_l`
  for (; p_diag_V_l < p_diag_V_l_max; ++p_scalarpot, ++p_rho_array, ++p_imagpot, ++p_diag_V_l) {
    rho = *p_rho_array;
    *p_diag_V_l = *p_scalarpot + coef / (rho * rho) - std::complex<double>(0.0,*p_imagpot);
  }
}

