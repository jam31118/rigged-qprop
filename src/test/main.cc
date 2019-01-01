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

//// headers for MPI
#include "mpi.h"

int error_and_exit(int rank, int return_code, const char *func_name) {
  fprintf(stderr, "In process with rank '%d': something got wrong in '%s' with return code: '%d'\n",
      rank, func_name, return_code);
  MPI_Finalize();
  return return_code;
}

//// `main` program starts here
int main(int argc, char *argv[]) {
  
  //// Declare variables with optional initialization
  int return_code = -1;
  int num_of_process = -1, rank = -1;
  
  //// MPI Initialization
  return_code = MPI_Init(&argc, &argv);
  if (return_code == MPI_SUCCESS) {
    fprintf(stdout, "[ LOG ] MPI program has been initialized.\n");
  } else {
    fprintf(stderr, "[ERROR] Something got wrong during initialization.\n");
    return -1;
  }

  //// Getting MPI communication size and process rank
  return_code = MPI_Comm_size(MPI_COMM_WORLD, &num_of_process);
  if (return_code != MPI_SUCCESS || num_of_process < 0) {
    fprintf(stderr, "[ERROR} Something got wrong during getting number of processes\n");
    MPI_Finalize();
    return -1;
  }
  return_code = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (return_code != MPI_SUCCESS || rank < 0) {
    fprintf(stderr, "[ERROR} Something got wrong during getting rank\n");
    MPI_Finalize();
    return -1;
  }
  
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
  // declare pointers for l-indenpendent, double arrays
  double 
    *rho_array = NULL, 
    *scalarpot_array = NULL, 
    *imagpot_array = NULL;
  // declare pointers for l-dependent, complex double arrays
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
  // open file
  string current_wf_bin_file_name = string("current-wf.bin");
  MPI_File fh;
  return_code = MPI_File_open(MPI_COMM_WORLD, current_wf_bin_file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_open"); }

  // Get file size
  MPI_Offset file_size;
  return_code = MPI_File_get_size(fh, &file_size);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_get_size"); }
  if (rank == 0) { fprintf(stdout, "[@rank=%d] The size of file '%s' = %lld\n", rank, current_wf_bin_file_name.c_str(), file_size); }
  long current_wf_file_size = file_size;

//  std::ifstream current_wf_bin_file;
//  string current_wf_bin_file_name = string("current-wf.bin");
//  current_wf_bin_file.open(current_wf_bin_file_name, std::ifstream::ate | std::ifstream::binary);
//  long current_wf_file_size = current_wf_bin_file.tellg();
//  current_wf_bin_file.close();
//  std::cout << "[ LOG ] current_wf_file_size: " << current_wf_file_size << std::endl;
  
  // Determine the number of basis wf (i.e. `wf_lm`) in the wf file
  long num_of_wf_lm = g_prop.num_of_phi_lm();
  if (num_of_wf_lm < 0) { std::cerr << "[ERROR] during `g_prop.num_of_phi_lm()`\n"; return -1; }
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
  long num_of_wf_per_proc;
  num_of_wf_per_proc = ( num_of_wf_lm + (num_of_process - 1) ) / num_of_process;
  if ( num_of_wf_lm >= num_of_process ) {
    num_of_wf_to_read = num_of_wf_per_proc;
    if (rank == num_of_process - 1) {
      num_of_wf_to_read = num_of_wf_lm - (num_of_process - 1) * num_of_wf_per_proc;
    } 
  } else {
      if (rank < num_of_wf_lm) { num_of_wf_to_read = 1; }
      else { num_of_wf_to_read = 0; }
  }
//  num_of_wf_to_read = num_of_wf_lm;
//  lm_index_start = 0;
  lm_index_start = rank * num_of_wf_per_proc;

  lm_index_max = lm_index_start + num_of_wf_to_read;
  if (lm_index_max > num_of_wf_lm) {
    std::cerr << "[ERROR] `lm_index_max` should not exceed `num_of_wf_lm`\n";
    return -1;
  }

  //// specify wf region to read
  long offset_lm;
//  MPI_Offset offset_lm;
  offset_lm = lm_index_start * bytes_per_wf;
  std::cout << "[ LOG ] offset_lm = " << offset_lm << std::endl;
  std::cout << "[ LOG ] bytes_per_wf = " << bytes_per_wf << std::endl;

  //// Read wf data from file
  std::complex<double> *wf_read = new std::complex<double>[N_rho * num_of_wf_to_read];
  MPI_Status read_status;
  return_code = MPI_File_read_at(fh, offset_lm, wf_read, num_of_wf_to_read * N_rho, MPI::DOUBLE_COMPLEX, &read_status);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_read_at"); }

  // Close file
  return_code = MPI_File_close(&fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_close"); }

//  current_wf_bin_file.open(current_wf_bin_file_name, std::ifstream::binary);
//  current_wf_bin_file.seekg(offset_lm);
//  long bytes_to_read = num_of_wf_to_read * bytes_per_wf;
//  current_wf_bin_file.read((char *)wf_read, bytes_to_read);
//  current_wf_bin_file.close();

  //// Allocate memory
  // declare pointers for l-indenpendent, double arrays
  rho_array = new double[N_rho];
  scalarpot_array = new double[N_rho];
  imagpot_array = new double[N_rho];
  // declare pointers for l-denpendent, double arrays
  diag_unitary = new std::complex<double>[N_rho];
  lower_offdiag_unitary = new std::complex<double>[N_rho-1];
  upper_offdiag_unitary = new std::complex<double>[N_rho-1];
  diag_unitary_inv = new std::complex<double>[N_rho];
  lower_offdiag_unitary_inv = new std::complex<double>[N_rho-1];
  upper_offdiag_unitary_inv = new std::complex<double>[N_rho-1];

  //// Instatiating potential objects
  // Prepare scalarpot
  scalarpot atomic_potential(
      para_ini.getDouble("nuclear-charge"), 
      para_ini.getDouble("pot-cutoff"), 
      get_effpot_alpha(para_ini));
  // Prepare imagpot
  const long imag_potential_width=long(para_prop.getDouble("imag-width")/delta_rho);
  imagpot imaginarypot(imag_potential_width);

  //// Prepare l-indenepent arrays
  double rho_value;
  for (i=0; i<N_rho; i++) {
    rho_value = delta_rho * (i+1);
    rho_array[i] = rho_value;
    scalarpot_array[i] = atomic_potential(rho_value, 0, 0, 0, 0);
    imagpot_array[i] = imaginarypot(i, 0, 0, 0, N_rho);
  }
  
  //// Allocate memory for storing propagators for each `l` values
  std::complex<double> *diag_unitary_stack = new std::complex<double>[N_rho * num_of_wf_to_read];
  std::complex<double> *lower_offdiag_unitary_stack = new std::complex<double>[(N_rho-1) * num_of_wf_to_read];
  std::complex<double> *upper_offdiag_unitary_stack = new std::complex<double>[(N_rho-1) * num_of_wf_to_read];
  std::complex<double> *diag_unitary_inv_stack = new std::complex<double>[N_rho * num_of_wf_to_read];
  std::complex<double> *lower_offdiag_unitary_inv_stack = new std::complex<double>[(N_rho-1) * num_of_wf_to_read];
  std::complex<double> *upper_offdiag_unitary_inv_stack = new std::complex<double>[(N_rho-1) * num_of_wf_to_read];

  //// Store propagators to each stacks
  long offset_in_diag_unitary_stack;
  long offset_in_offidag_unitary_stack;
  long diag_length;
  long offdiag_length;
  diag_length = N_rho;
  offdiag_length = N_rho-1;
  //// Start storing
  for (lm_index=lm_index_start; lm_index<lm_index_max; ++lm_index) {

    //// Retrieve `l` and `m` quantum numbers
    if (0 != get_ell_and_m_from_lm_index(lm_index, &l, &m, para_ini.getLong("initial-m"), g_prop.dimens())) { 
      fprintf(stderr, "[ERROR] during `get_ell_and_m_from_lm_index`\n");
      return 1; 
    }
    
    //// Evaulate unitary propagator
    // sign = -1 for explicit half time propagation
    sign = -1;
    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(
        l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,
        diag_unitary,lower_offdiag_unitary,upper_offdiag_unitary);
    // sign = 1 for implicit half time propagation
    sign = 1;
    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(
        l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,
        diag_unitary_inv,lower_offdiag_unitary_inv,upper_offdiag_unitary_inv);

    //// Store it to stack
    offset_in_diag_unitary_stack = (lm_index-lm_index_start)*N_rho;
    offset_in_offidag_unitary_stack = (lm_index-lm_index_start)*(N_rho-1);
    std::copy(diag_unitary, diag_unitary+diag_length, 
        diag_unitary_stack+offset_in_diag_unitary_stack);
    std::copy(lower_offdiag_unitary, lower_offdiag_unitary+offdiag_length, 
        lower_offdiag_unitary_stack+offset_in_offidag_unitary_stack);
    std::copy(upper_offdiag_unitary, upper_offdiag_unitary+offdiag_length, 
        upper_offdiag_unitary_stack+offset_in_offidag_unitary_stack);
    std::copy(diag_unitary_inv, diag_unitary_inv+diag_length, 
        diag_unitary_inv_stack+offset_in_diag_unitary_stack);
    std::copy(lower_offdiag_unitary_inv, lower_offdiag_unitary_inv+offdiag_length, 
        lower_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack);
    std::copy(upper_offdiag_unitary_inv, upper_offdiag_unitary_inv+offdiag_length, 
        upper_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack);
  } 

  std::cout << "[ LOG ] all unitary propagator stored\n";

  //// Set time index range for propagation
  long start_time_index = 0;
  double post_prop_duration = para_tsurff.getDouble("R-tsurff") / para_tsurff.getDouble("p-min-tsurff");
  double p_min_tsurff;
  long num_of_time_steps;
//  num_of_time_steps = long(post_prop_duration / delta_t);
  try { num_of_time_steps = para_prop.getLong("num-of-time-steps"); }
  catch (std::exception&) { 
    try { p_min_tsurff = para_tsurff.getDouble("p-min-tsurff"); }
    catch (std::exception&) {
      std::cerr << "[ERROR] either `num-of-time-steps` or `p-min-tsurff` should be set\n";
      return -1;
    }
    post_prop_duration = para_tsurff.getDouble("R-tsurff") / p_min_tsurff;
    num_of_time_steps = long(post_prop_duration / delta_t); 
  }
//  num_of_time_steps = 8000;

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

    //// Set offset in propagators stacks to select for the current `l` values
    offset_in_diag_unitary_stack = (lm_index-lm_index_start)*N_rho;
    offset_in_offidag_unitary_stack = (lm_index-lm_index_start)*(N_rho-1);
  
    //// Propagate
    // explicit half time propagation
    mat_vec_mul_tridiag(
        diag_unitary_stack+offset_in_diag_unitary_stack,
        lower_offdiag_unitary_stack+offset_in_offidag_unitary_stack,
        upper_offdiag_unitary_stack+offset_in_offidag_unitary_stack,
        wf_lm, wf_lm_mid, N_rho);
    // implicit half time propagation
    gaussian_elimination_tridiagonal(
        diag_unitary_inv_stack+offset_in_diag_unitary_stack,
        lower_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack,
        upper_offdiag_unitary_inv_stack+offset_in_offidag_unitary_stack,
        wf_lm, wf_lm_mid, N_rho);

    //// Evaulate tsurff-related values
    //
    //// end
      
    }
    //// Logging
    if (((time_index + 1) % num_of_steps_to_print_progress) == 0) {
      num_of_steps_done_so_far = time_index - start_time_index;
      if (rank == 0) {
        std::cout << "[@rank=" << rank << "][ LOG ] num_of_steps_done_so_far = " << num_of_steps_done_so_far + 1 << " / " << num_of_time_steps << std::endl;
      }
    }
  }
  std::cout << "[ LOG ] Propagation done\n";

  //// Write wf data to a file
//  MPI_File fh;
  return_code = MPI_File_open(MPI_COMM_WORLD, current_wf_bin_file_name.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_open"); }

  MPI_Status write_status;
  return_code = MPI_File_write_at(fh, offset_lm, wf_read, num_of_wf_to_read * N_rho, MPI::DOUBLE_COMPLEX, &write_status);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_write_at"); }

  //// Check `write_status` to confirm how many elements has been written

  // Close file
  return_code = MPI_File_close(&fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_close"); }

//  std::ofstream current_wf_bin_file_out;
//  current_wf_bin_file_out.open(current_wf_bin_file_name, std::ios::binary);
//  current_wf_bin_file_out.write((char *) wf_read, bytes_to_read);
//  current_wf_bin_file_out.close();
  
  //// Write propagator arrays to files
//  if (rank==0) {
//    FILE *diag_unitary_file = fopen("diag_unitary_stack.bin", "wb");
//    fwrite(diag_unitary_stack, sizeof(std::complex<double>), N_rho*num_of_wf_lm, diag_unitary_file);
//    fclose(diag_unitary_file);
//    FILE *lower_offdiag_unitary_file = fopen("lower_offdiag_unitary_stack.bin", "wb");
//    fwrite(lower_offdiag_unitary_stack, sizeof(std::complex<double>), (N_rho-1)*num_of_wf_lm, lower_offdiag_unitary_file);
//    fclose(lower_offdiag_unitary_file);
//    FILE *upper_offdiag_unitary_file = fopen("upper_offdiag_unitary_stack.bin", "wb");
//    fwrite(upper_offdiag_unitary_stack, sizeof(std::complex<double>), (N_rho-1)*num_of_wf_lm, upper_offdiag_unitary_file);
//    fclose(upper_offdiag_unitary_file);
//  }
  
  //// MPI Finalization
  return_code = MPI_Finalize();
  if (return_code == MPI_SUCCESS) {
    fprintf(stdout, "[ LOG ] MPI program has been finalized.\n");
  } else {
    fprintf(stderr, "[ERROR] Something got wrong during finalization.\n");
    return -1;
  }

  //// End this program
  return 0;
}
