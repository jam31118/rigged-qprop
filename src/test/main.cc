//// headers from standard libraries
#include <iostream>
#include <complex>
#include <string>
#include <fstream>
#include <algorithm>

//// headers from `QPROP`
#include "potentials.hh"
#include "parameter.hh"
#include "common.hh"
#include "propagator.hh"

//// headers from `matrix`
#include "matrix.hh"

//// `main` program starts here
int main() {

  //// Load parameter files for QPROP
  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");

  //// Declare variables
  grid g_prop;
  long i; // for iteration
  long sign; // a sign in unitary propagator (-1: forward, 1: backward propagation)
  long l, m, N_rho;
  long lm_index; // iteration for lm_unifired index
  double Z, delta_rho, delta_t, grid_size;
  double 
    *rho_array = NULL, 
    *scalarpot_array = NULL, 
    *imagpot_array = NULL;
  std::complex<double> 
    *diag_unitary = NULL, 
    *lower_offdiag_unitary = NULL, 
    *upper_offdiag_unitary = NULL,
    *diag_unitary_inv = NULL, 
    *lower_offdiag_unitary_inv = NULL, 
    *upper_offdiag_unitary_inv = NULL;
  
  //// Extract calculation parameters
  delta_rho = para_ini.getDouble("delta-r");
  grid_size = get_grid_size(para_ini, para_prop, para_tsurff);
  N_rho = long(grid_size/delta_rho);
  Z = para_ini.getDouble("nuclear-charge");
  delta_t = para_prop.getDouble("delta-t");
  
  //// Set grid parameters
  g_prop.set_dim(para_ini.getLong("qprop-dim"));
  g_prop.set_ngps(N_rho, para_prop.getLong("ell-grid-size"), 1);
  g_prop.set_delt(delta_rho);
  g_prop.set_offs(0, 0, 0);

  //// Configure wf file 
  std::ifstream current_wf_bin_file;
  string current_wf_bin_file_name = string("current-wf.bin");
  current_wf_bin_file.open(current_wf_bin_file_name, std::ifstream::ate | std::ifstream::binary);
  long current_wf_file_size = current_wf_bin_file.tellg();
  current_wf_bin_file.close();
  std::cout << "[ LOG ] current_wf_file_size: " << current_wf_file_size << std::endl;
  
  // Dertermine the number of basis wf (i.e. `wf_lm`) in the wf file
  long num_of_wf_lm;
  switch (g_prop.dimens()) {
    case 34: num_of_wf_lm = g_prop.ngps_y(); break;
    case 44: num_of_wf_lm = g_prop.ngps_y() * g_prop.ngps_y(); break;
    default:
             std::cerr << "[ERROR] unknown qprop-dimension: " 
               << g_prop.dimens() << std::endl;
             return -1;
  }
  std::cout << "[ LOG ] num_of_wf_lm = " << num_of_wf_lm << std::endl;
  // Determine the `bytes_per_wf` to read
  long bytes_per_wf = current_wf_file_size / num_of_wf_lm;
  // Determine the number of numbers per wf and check if it is same as `N_rho`
  long num_of_numbers_per_wf = bytes_per_wf / sizeof(std::complex<double>);
  if (num_of_numbers_per_wf != N_rho) {
    fprintf(stderr, "[ERROR] num_of_numbers_per_wf != N_rho\n");
    return -1;
  }

  //// Set parameter
  long lm_index_start, lm_index_max, num_of_wf_to_read;
  lm_index_start = 0;
  num_of_wf_to_read = num_of_wf_lm;
  lm_index_max = lm_index_start + num_of_wf_to_read;
  if (lm_index_max > num_of_wf_lm) {
    std::cerr << "[ERROR] `lm_index_max` should not exceed `num_of_wf_lm`\n";
    return -1;
  }

  //// specify wf region to read
  long offset_lm;
  offset_lm = lm_index_start * bytes_per_wf;
  std::cout << "[ LOG ] offset_lm = " << offset_lm << std::endl;
  std::cout << "[ LOG ] bytes_per_wf = " << bytes_per_wf << std::endl;

  //// Read wf data from file
  std::complex<double> *wf_read = new std::complex<double>[N_rho * num_of_wf_to_read];
  current_wf_bin_file.open(current_wf_bin_file_name, std::ifstream::binary);
  current_wf_bin_file.seekg(offset_lm);
  long bytes_to_read = num_of_wf_to_read * bytes_per_wf;
  current_wf_bin_file.read((char *)wf_read, bytes_to_read);
  current_wf_bin_file.close();

 
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
  
  //// Allocate memory for propagator
  std::complex<double> *diag_unitary_stack = new std::complex<double>[N_rho * num_of_wf_to_read];
  std::complex<double> *lower_offdiag_unitary_stack = new std::complex<double>[(N_rho-1) * num_of_wf_to_read];
  std::complex<double> *upper_offdiag_unitary_stack = new std::complex<double>[(N_rho-1) * num_of_wf_to_read];
  std::complex<double> *diag_unitary_inv_stack = new std::complex<double>[N_rho * num_of_wf_to_read];
  std::complex<double> *lower_offdiag_unitary_inv_stack = new std::complex<double>[(N_rho-1) * num_of_wf_to_read];
  std::complex<double> *upper_offdiag_unitary_inv_stack = new std::complex<double>[(N_rho-1) * num_of_wf_to_read];

  long offset_in_diag_unitary_stack;
  long offset_in_offidag_unitary_stack;
  long diag_length;
  long offdiag_length;
//  std::complex<double> *unitary_prop_lm, *p_unitary_prop_lm;
  for (lm_index=lm_index_start; lm_index<lm_index_max; ++lm_index) {

    //// Retrieve `l` and `m` quantum numbers
    if (0 != get_ell_and_m_from_lm_index(lm_index, &l, &m, para_ini.getLong("initial-m"), g_prop.dimens())) { 
      fprintf(stderr, "[ERROR] during `get_ell_and_m_from_lm_index`\n");
      return 1; 
    }
    
    //// Evaulate unitary propagator
    sign = -1;
    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,diag_unitary,lower_offdiag_unitary,upper_offdiag_unitary);
    sign = 1;
    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,diag_unitary_inv,lower_offdiag_unitary_inv,upper_offdiag_unitary_inv);

    //// Store it to stack
    offset_in_diag_unitary_stack = (lm_index-lm_index_start)*N_rho;
    offset_in_offidag_unitary_stack = (lm_index-lm_index_start)*(N_rho-1);
    diag_length = N_rho;
    offdiag_length = N_rho-1;

    std::cout << "[ LOG ] before `copy()` at lm_index: " << lm_index << std::endl;

    std::copy(diag_unitary, diag_unitary+diag_length, diag_unitary_stack+offset_in_diag_unitary_stack);
    std::copy(lower_offdiag_unitary, lower_offdiag_unitary+offdiag_length, lower_offdiag_unitary_stack+offset_in_offidag_unitary_stack);
    std::copy(upper_offdiag_unitary, upper_offdiag_unitary+offdiag_length, upper_offdiag_unitary_stack+offset_in_offidag_unitary_stack);
    std::copy(diag_unitary_inv, diag_unitary_inv+diag_length, diag_unitary_inv_stack+offset_in_diag_unitary_stack);
    std::copy(lower_offdiag_unitary_inv, lower_offdiag_unitary_inv+offdiag_length, lower_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack);
    std::copy(upper_offdiag_unitary_inv, upper_offdiag_unitary_inv+offdiag_length, upper_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack);
//    for (i=offset; i<offset+N_rho; ++i) { diag_unitary_stack[i] = }
//    for (p_unitary_prop_lm=unitary_prop_lm; p_unitary_prop_lm<p_unitary_prop_lm_max; ++p_unitary_prop_lm) {
//      *p_unitary_prop_lm = ;
//    }
  } 

  std::cout << "[ LOG ] all unitary propagator stored\n";

  //// Set time index range for propagation
  long start_time_index = 0;
  double post_prop_duration = para_tsurff.getDouble("R-tsurff") / para_tsurff.getDouble("p-min-tsurff");
  long num_of_time_steps;
//  num_of_time_steps = long(post_prop_duration / delta_t);
  num_of_time_steps = 8000;

  long time_index_max = start_time_index + num_of_time_steps;
  long time_index;
  long num_of_steps_to_print_progress = 200; // [NOTE] to become global default config
  std::cout << "[ LOG ] num_of_time_steps: " << num_of_time_steps << std::endl;
  std::cout << "[ LOG ] post_prop_duration: " << post_prop_duration << std::endl;

  //// Start iteration over each `wf_lm` for `lm_index`
  long num_of_steps_done_so_far;
  long num_of_numbers_before_this_wf_lm;
  std::complex<double> *wf_lm = wf_read;
  std::complex<double> *wf_lm_mid = new std::complex<double>[N_rho];
  for (time_index=start_time_index; time_index<time_index_max; ++time_index) {
    for (lm_index=lm_index_start; lm_index<lm_index_max; ++lm_index) {

    //// Set wf_lm pointer
    num_of_numbers_before_this_wf_lm = (lm_index - lm_index_start) * num_of_numbers_per_wf;
    wf_lm = wf_read + num_of_numbers_before_this_wf_lm; 
    
//    //// Retrieve `l` and `m` quantum numbers
//    if (0 != get_ell_and_m_from_lm_index(lm_index, &l, &m, g_prop.dimens())) { 
//      fprintf(stderr, "[ERROR] during `get_ell_and_m_from_lm_index`\n");
//      return 1; 
//    };

    offset_in_diag_unitary_stack = (lm_index-lm_index_start)*N_rho;
    offset_in_offidag_unitary_stack = (lm_index-lm_index_start)*(N_rho-1);
  
    //// Evaulate unitary propagator
//    sign = -1;
//    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,diag_unitary_stack+offset_in_diag_unitary_stack,lower_offdiag_unitary_stack+offset_in_offidag_unitary_stack,upper_offdiag_unitary_stack+offset_in_offidag_unitary_stack);
//    sign = 1;
//    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,diag_unitary_inv_stack+offset_in_diag_unitary_stack,lower_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack,upper_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack);
  
    //// Propagate
    mat_vec_mul_tridiag(diag_unitary_stack+offset_in_diag_unitary_stack,lower_offdiag_unitary_stack+offset_in_offidag_unitary_stack,upper_offdiag_unitary_stack+offset_in_offidag_unitary_stack, wf_lm, wf_lm_mid, N_rho);
    gaussian_elimination_tridiagonal(diag_unitary_inv_stack+offset_in_diag_unitary_stack,lower_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack,upper_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack, wf_lm, wf_lm_mid, N_rho);
    //// Evaulate tsurff-related values
    //
    //// end
      
    }
    //// Logging
    if (((time_index + 1) % num_of_steps_to_print_progress) == 0) {
      num_of_steps_done_so_far = time_index - start_time_index;
      std::cout << "[ LOG ] num_of_steps_done_so_far = " << num_of_steps_done_so_far + 1 << " / " << num_of_time_steps << std::endl;
//      std::cout << "[ LOG ] time_index = " << time_index << " | " << "((time_index + 1) % num_of_steps_to_print_progress) == 0: " << (((time_index + 1) % num_of_steps_to_print_progress) == 0) << std::endl;
//      std::cout << "[ LOG ] wf_lm : wf_lm + N_rho : wf_read_max = " << wf_lm << " : " << (wf_lm + N_rho) << " : " << (wf_read + num_of_numbers_per_wf * num_of_wf_to_read) << std::endl;
    }
  }
  std::cout << "[ LOG ] Propagation done\n";


  //// Write wf data to file
  std::ofstream current_wf_bin_file_out;
  current_wf_bin_file_out.open(current_wf_bin_file_name, std::ios::binary);
  current_wf_bin_file_out.write((char *) wf_read, bytes_to_read);
  current_wf_bin_file_out.close();
  

  //// Print results
  FILE *diag_unitary_file = fopen("diag_unitary_stack.bin", "wb");
  fwrite(diag_unitary_stack, sizeof(std::complex<double>), N_rho*num_of_wf_lm, diag_unitary_file);
  fclose(diag_unitary_file);
  FILE *lower_offdiag_unitary_file = fopen("lower_offdiag_unitary_stack.bin", "wb");
  fwrite(lower_offdiag_unitary_stack, sizeof(std::complex<double>), (N_rho-1)*num_of_wf_lm, lower_offdiag_unitary_file);
  fclose(lower_offdiag_unitary_file);
  FILE *upper_offdiag_unitary_file = fopen("upper_offdiag_unitary_stack.bin", "wb");
  fwrite(upper_offdiag_unitary_stack, sizeof(std::complex<double>), (N_rho-1)*num_of_wf_lm, upper_offdiag_unitary_file);
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
