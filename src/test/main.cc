#include <iostream>
#include <complex>
#include <string>
#include <fstream>
#include "potentials.hh"
#include "parameter.hh"
#include "common.hh"
#include "propagator.hh"

#include "matrix.hh"

int main() {

  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");

  //// Declare variables
  
  grid g_prop;

  long i; // for iteration
  long l, sign;
  long N_rho;
  double Z, delta_rho, delta_t;
  double *rho_array = NULL, *scalarpot_array = NULL, *imagpot_array = NULL;
  std::complex<double> 
    *diag_unitary = NULL, 
    *lower_offdiag_unitary = NULL, 
    *upper_offdiag_unitary = NULL,
    *diag_unitary_inv = NULL, 
    *lower_offdiag_unitary_inv = NULL, 
    *upper_offdiag_unitary_inv = NULL;
  
  //// Extract calculation parameters
  delta_rho = para_ini.getDouble("delta-r");
  double grid_size = get_grid_size(para_ini, para_prop, para_tsurff);
  N_rho = long(grid_size/delta_rho);
  Z = para_ini.getDouble("nuclear-charge");
  delta_t = para_prop.getDouble("delta-t");
  
  //// Set grid parameters
  g_prop.set_dim(para_ini.getLong("qprop-dim"));
  g_prop.set_ngps(N_rho, para_prop.getLong("ell-grid-size"), 1);
  g_prop.set_delt(delta_rho);
  g_prop.set_offs(0, 0, 0);

  //// Set parameter
  long lm_index = 2;
  long m;
  if (0 != get_ell_and_m_from_lm_index(lm_index, &l, &m)) { 
    fprintf(stderr, "[ERROR] during `get_ell_and_m_from_lm_index`\n");
    return 1; 
  };
//  l = 1;

  //// Configure wf file 
  std::ifstream current_wf_bin_file;
  string current_wf_bin_file_name = string("current-wf.bin");
  current_wf_bin_file.open(current_wf_bin_file_name, std::ifstream::ate | std::ifstream::binary);
  long current_wf_file_size = current_wf_bin_file.tellg();
  current_wf_bin_file.close();
  std::cout << "[ LOG ] current_wf_file_size: " << current_wf_file_size << std::endl;

  //// specify wf region to read
  long num_of_wf_lm;
  if (g_prop.dimens() == 34) {num_of_wf_lm = g_prop.ngps_y();}
  else if (g_prop.dimens() == 44) {num_of_wf_lm = g_prop.ngps_y() * g_prop.ngps_y();}

  std::cout << "[ LOG ] num_of_wf_lm = " << num_of_wf_lm << std::endl;
  long bytes_per_wf = current_wf_file_size / num_of_wf_lm;
  long num_of_numbers_per_wf = bytes_per_wf / sizeof(std::complex<double>);
  if (num_of_numbers_per_wf != N_rho) {
    fprintf(stderr, "[ERROR] num_of_numbers_per_wf != N_rho\n");
    return -1;
  }
  long offset_lm = lm_index * bytes_per_wf;
  std::cout << "[ LOG ] offset_lm = " << offset_lm << std::endl;
  std::cout << "[ LOG ] bytes_per_wf = " << bytes_per_wf << std::endl;

  //// Read wf data from file
  std::complex<double> *wf_lm = new std::complex<double>[N_rho];
  current_wf_bin_file.open(current_wf_bin_file_name, std::ifstream::binary);
  current_wf_bin_file.seekg(offset_lm);
  current_wf_bin_file.read((char *)wf_lm, bytes_per_wf);
  current_wf_bin_file.close();
 
//  for (i=0; i<N_rho; i++) {
//    std::cout << wf_lm[i] << " ";
//  } std::cout << std::endl;

  //// Allocate memory
  rho_array = new double[N_rho];
  scalarpot_array = new double[N_rho];
  imagpot_array = new double[N_rho];

  diag_unitary = new std::complex<double>[N_rho];
  lower_offdiag_unitary = new std::complex<double>[N_rho-1];
  upper_offdiag_unitary = new std::complex<double>[N_rho-1];
  
  diag_unitary_inv = new std::complex<double>[N_rho];
  lower_offdiag_unitary_inv = new std::complex<double>[N_rho-1];
  upper_offdiag_unitary_inv = new std::complex<double>[N_rho-1];
//  w_array = new std::complex<double>[N_rho];

  //// Instatiating potential objects
  scalarpot atomic_potential(para_ini.getDouble("nuclear-charge"), para_ini.getDouble("pot-cutoff"), get_effpot_alpha(para_ini));
  const long imag_potential_width=long(para_prop.getDouble("imag-width")/delta_rho);
  imagpot imaginarypot(imag_potential_width);

  //// Prepare arrays
  double rho_value;
  for (i=0; i<N_rho; i++) {
    rho_value = delta_rho * (i+1);
    rho_array[i] = rho_value;
    scalarpot_array[i] = atomic_potential(rho_value, 0, 0, 0, 0);
    imagpot_array[i] = imaginarypot(i, 0, 0, 0, N_rho);
  }

  //// Evaulate
  sign = -1;
  evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,diag_unitary,lower_offdiag_unitary,upper_offdiag_unitary);
  sign = 1;
  evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,diag_unitary_inv,lower_offdiag_unitary_inv,upper_offdiag_unitary_inv);

  //// Propagate
  std::complex<double> *wf_lm_mid = new std::complex<double>[N_rho];
  long start_time_index = 0;
  long num_of_time_steps = 5;
  long time_index_max = start_time_index + num_of_time_steps;
  long time_index;
  for (time_index=start_time_index; time_index<time_index_max; time_index++) {
    mat_vec_mul_tridiag(diag_unitary, lower_offdiag_unitary, upper_offdiag_unitary, wf_lm, wf_lm_mid, N_rho);
    gaussian_elimination_tridiagonal(diag_unitary_inv, lower_offdiag_unitary_inv, upper_offdiag_unitary_inv, wf_lm, wf_lm_mid, N_rho);
    //// Evaulate tsurff-related values
    //
    //// end
  }

  std::cout << "[ LOG ] Propgation done\n";

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


//  for (i=0; i<N_rho; i++) {
//    std::cout << diag_unitary[i] << " ";
//  } std::cout << std::endl;
//  for (i=0; i<N_rho-1; i++) {
//    std::cout << lower_offdiag_unitary[i] << " ";
//  } std::cout << std::endl;
//  for (i=0; i<N_rho-1; i++) {
//    std::cout << upper_offdiag_unitary[i] << " ";
//  } std::cout << std::endl;
  
  //// End this program
  return 0;
}
