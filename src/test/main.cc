#include <iostream>
#include <complex>
#include "potentials.hh"
#include "parameter.hh"
#include "common.hh"
#include "propagator.hh"

int main() {

  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");

  //// Declare variables
  long l, sign;
  long N_rho;
  double Z, delta_rho, delta_t;
  double *rho_array = NULL, *scalarpot_array = NULL, *imagpot_array = NULL;
  std::complex<double> *diag_unitary = NULL, *lower_offdiag_unitary = NULL, *upper_offdiag_unitary = NULL; // *w_array = NULL;
  
  //// Extract calculation parameters
  delta_rho = para_ini.getDouble("delta-r");
  double grid_size = get_grid_size(para_ini, para_prop, para_tsurff);
  N_rho = long(grid_size/delta_rho);
  Z = para_ini.getDouble("nuclear-charge");
  delta_t = para_prop.getDouble("delta-t");

  //// Allocate memory
  rho_array = new double[N_rho];
  scalarpot_array = new double[N_rho];
  imagpot_array = new double[N_rho];
  diag_unitary = new std::complex<double>[N_rho];
  lower_offdiag_unitary = new std::complex<double>[N_rho-1];
  upper_offdiag_unitary = new std::complex<double>[N_rho-1];
//  w_array = new std::complex<double>[N_rho];

  //// Instatiating potential objects
  scalarpot atomic_potential(para_ini.getDouble("nuclear-charge"), para_ini.getDouble("pot-cutoff"), get_effpot_alpha(para_ini));
  const long imag_potential_width=long(para_prop.getDouble("imag-width")/delta_rho);
  imagpot imaginarypot(imag_potential_width);

  //// Prepare arrays
  double rho_value;
  long i;
  for (i=0; i<N_rho; i++) {
    rho_value = delta_rho * (i+1);
    rho_array[i] = rho_value;
    scalarpot_array[i] = atomic_potential(rho_value, 0, 0, 0, 0);
    imagpot_array[i] = imaginarypot(i, 0, 0, 0, N_rho);
  }

  //// Set parameter
  l = 1;
  sign = -1;

  //// Evaulate
  evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,diag_unitary,lower_offdiag_unitary,upper_offdiag_unitary);

  //// Print results
  FILE *diag_unitary_file = fopen("diag_unitary.bin", "wb");
  fwrite(diag_unitary, sizeof(std::complex<double>), N_rho, diag_unitary_file);
  fclose(diag_unitary_file);
  FILE *lower_offdiag_unitary_file = fopen("lower_offdiag_unitary.bin", "wb");
  fwrite(lower_offdiag_unitary, sizeof(std::complex<double>), N_rho-1, lower_offdiag_unitary_file);
  fclose(lower_offdiag_unitary_file);
  FILE *upper_offdiag_unitary_file = fopen("upper_offdiag_unitary.bin", "wb");
  fwrite(upper_offdiag_unitary, sizeof(std::complex<double>), N_rho-1, upper_offdiag_unitary_file);
  fclose(upper_offdiag_unitary_file);


  for (i=0; i<N_rho; i++) {
    std::cout << diag_unitary[i] << " ";
  } std::cout << std::endl;
  for (i=0; i<N_rho-1; i++) {
    std::cout << lower_offdiag_unitary[i] << " ";
  } std::cout << std::endl;
  for (i=0; i<N_rho-1; i++) {
    std::cout << upper_offdiag_unitary[i] << " ";
  } std::cout << std::endl;
  
  //// End this program
  return 0;
}
