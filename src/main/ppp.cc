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
#include "tridiag-common.hh"

//// header(s) for MPI
#ifdef HAVE_MPI 
#include "mpi.h"
#endif

//// header(s) for `ppp-gpu`
#ifdef HAVE_CUDA
#include "cu_propagator.h"
#endif

//// header for time-measurement
#include <chrono>
//#ifdef HAVE_BOOST
//#include <boost/timer.hpp>
//#endif

//// For time measurement
//typedef std::chrono::high_resolution_clock Clock;

// Define macro function
#define MIN(x,y) ((x<y)?x:y)
//#define DUMMY_INT -1<<31


int error_and_exit(int rank, int return_code, const char *func_name) {
  fprintf(stderr, "In process with rank '%d': something got wrong in '%s' with return code: '%d'\n",
      rank, func_name, return_code);
#ifdef HAVE_MPI 
  MPI_Finalize();
#endif
  return return_code;
}


#ifdef HAVE_MPI 
int finalize_mpi_with_check(int rank) {
  int return_code = MPI_Finalize();
  if (return_code == MPI_SUCCESS) {
    fprintf(stdout, "[@rank=%d][ LOG ] MPI program has been finalized.\n", rank);
  } else {
    fprintf(stderr, "[ERROR] Something got wrong during finalization.\n");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
#endif


//// `main` program starts here
int main(int argc, char *argv[]) {
  
  //// Declare variables with optional initialization
  int return_code = EXIT_FAILURE;
  int num_of_process = 1, rank = 0;  // default is a single-process mode
  
#ifdef HAVE_MPI 
  //// MPI Initialization
  return_code = MPI_Init(&argc, &argv);
  if (return_code == MPI_SUCCESS) {
    fprintf(stdout, "[ LOG ] MPI program has been initialized.\n");
  } else {
    fprintf(stderr, "[ERROR] Something got wrong during initialization.\n");
    return return_code;
  }
#endif


#ifdef HAVE_MPI 
  //// Getting MPI communication size and process rank
  return_code = MPI_Comm_size(MPI_COMM_WORLD, &num_of_process);
  if (return_code != MPI_SUCCESS || num_of_process < 0) {
    fprintf(stderr, "[ERROR} Something got wrong during getting number of processes\n");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  return_code = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (return_code != MPI_SUCCESS || rank < 0) {
    fprintf(stderr, "[ERROR} Something got wrong during getting rank\n");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
#endif


#ifdef HAVE_MPI

  //// MPI general configrations

  // Set `element_type`
  MPI_Datatype element_type;
  element_type = MPI::DOUBLE_COMPLEX;

  // Get `element_type`s size
  int element_type_size;
  MPI_Type_size(element_type, &element_type_size);

#endif



  //// Load parameter files for QPROP
  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");
  


  //// Get and set PPP configuation
  bool use_ppp = false;
  try { use_ppp = para_tsurff.getBool("use-ppp"); }
  catch (std::exception&) {}
  // Abort if `ppp` is not configured to run
  if (!use_ppp) { 
    fprintf(stderr, "[ERROR][@rank=%d] `ppp` hasn't been configured to run.\n", rank);
    fprintf(stderr, "[ERROR][@rank=%d] Please set `use-ppp` in `tsurff.param` file to `1` to run it.\n", rank);
    return EXIT_FAILURE; 
  }



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



  // Determine the number of basis wf (i.e. `wf_lm`) in the wf file
  long num_of_wf_lm = g_prop.num_of_phi_lm();
  if (num_of_wf_lm < 0) { std::cerr << "[ERROR] during `g_prop.num_of_phi_lm()`\n"; return EXIT_FAILURE; }
  if (rank == 0) {
    fprintf(stdout, "[ LOG ][@rank=%d] num_of_wf_lm = %d\n", rank, num_of_wf_lm);
  }
//  std::cout << "[ LOG ] num_of_wf_lm = " << num_of_wf_lm << std::endl;
  const int num_of_remaining_wf = num_of_wf_lm % num_of_process;



  //// Drop process that may be unused
  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group working_group; // = MPI_GROUP_EMPTY;
  int working_rank_range_incl[1][3] = { { 0, num_of_process-1, 1 } };
//  int working_rank_range_excl[1][3] = { { 0, num_of_process-1, 1 } };
  if ( num_of_wf_lm < num_of_process ) {
//    working_rank_range_excl[0][0] = num_of_remaining_wf;
    working_rank_range_incl[0][1] = num_of_remaining_wf-1;
//    MPI_Group_range_excl(working_group, 1, working_rank_range_excl, &working_group);
  }
  MPI_Group_range_incl(world_group, 1, working_rank_range_incl, &working_group);

  MPI_Comm working_world;
  MPI_Comm_create(MPI_COMM_WORLD, working_group, &working_world);

  if (working_world == MPI_COMM_NULL) {
    fprintf(stdout, "[ LOG ][@rank=%d] I'm going to finalize MPI\n", rank);
    MPI_Finalize();
    return EXIT_SUCCESS;
  }





  //// Configure wf file 
  // open file
  string current_wf_bin_file_name = string("current-wf.bin");
#ifdef HAVE_MPI
  MPI_File fh;
  return_code = MPI_File_open(working_world, current_wf_bin_file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_open"); }
#else
  std::fstream fh;
  fh.open(current_wf_bin_file_name, std::ios::in | std::ios::binary);
  if (!fh.is_open()) { fprintf(stderr,"[ERROR] during opening file `%s`\n", current_wf_bin_file_name.c_str()); }
#endif


  // Get file size
#ifdef HAVE_MPI
  MPI_Offset file_size;
  return_code = MPI_File_get_size(fh, &file_size);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_get_size"); }
  if (rank == 0) { fprintf(stdout, "[@rank=%d] The size of file '%s' = %lld\n", rank, current_wf_bin_file_name.c_str(), file_size); }
#else
  long file_size;
  fh.seekg(0, fh.end);
  file_size = fh.tellg();
  fh.seekg(0, fh.beg);
#endif
  long current_wf_file_size = file_size;


  // Determine the `bytes_per_wf` to read
  long bytes_per_wf = current_wf_file_size / num_of_wf_lm;
  if (rank == 0) {
    std::cout << "[@rank=" << rank << "]" << "[ LOG ] bytes_per_wf = " << bytes_per_wf << std::endl;
  }

  // Determine the number of numbers per wf and check if it is same as `N_rho`
  long num_of_numbers_per_wf = bytes_per_wf / sizeof(std::complex<double>);
  if (num_of_numbers_per_wf != N_rho) {
    fprintf(stderr, "[ERROR] num_of_numbers_per_wf (=%ld) != N_rho (=%ld) \n", num_of_numbers_per_wf, N_rho);
    return EXIT_FAILURE;
  }


  //// Set parameter
//  int lm_index_start, lm_index_max;
//  int num_of_wf_to_read;
  const int num_of_wf_per_proc_max = ( num_of_wf_lm + (num_of_process - 1) ) / num_of_process;
  const int num_of_wf_per_proc_min = num_of_wf_lm / num_of_process;
//  if ( num_of_wf_lm >= num_of_process ) {
//    num_of_wf_to_read = num_of_wf_per_proc;
//    if (rank == num_of_process - 1) {
//      num_of_wf_to_read = num_of_wf_lm - (num_of_process - 1) * num_of_wf_per_proc;
//    } 
//  } else {
//      if (rank < num_of_wf_lm) { num_of_wf_to_read = 1; }
//      else { num_of_wf_to_read = 0; }
//  }

  const int num_of_wf_to_read = num_of_wf_per_proc_min + (rank < num_of_remaining_wf);
//  const int num_of_wf_before_this_rank = ( (num_of_wf_lm / num_of_process) + (num_of_remaining_wf != 0) ) * rank - (rank - num_of_remaining_wf) * (rank >= num_of_remaining_wf);

//  lm_index_start = rank * num_of_wf_per_proc;
//  lm_index_start = num_of_wf_before_this_rank;
  int lm_index_start = rank;
//  if (rank >= num_of_remaining_wf) { lm_index_start = DUMMY_INT; }
//  const int lm_index_max = lm_index_start + num_of_wf_to_read;
  int lm_index_max = lm_index_start + num_of_wf_to_read * num_of_process;
//  if (rank >= num_of_remaining_wf) { lm_index_max = DUMMY_INT; }

  if ( (num_of_wf_lm < num_of_process) && (rank >= num_of_remaining_wf) ) {
    if ( (num_of_wf_to_read != 0) || (lm_index_max != lm_index_start) ) {
      fprintf(stderr, "[ERROR][@rank=%d] There is an inconsistency\n", rank);
      return EXIT_FAILURE;
    }
//    fprintf(stdout, "[ LOG ][@rank=%d] going to finalize MPI\n", rank);
//    MPI_Finalize();
//    return EXIT_SUCCESS;
  }
//  if ( (num_of_wf_lm < num_of_process) && (rank >= num_of_remaining_wf) ) {
//    fprintf(stdout, "[ LOG ][@rank=%d] going to finalize MPI\n", rank);
//    MPI_Finalize();
//    fprintf(stdout, "[ LOG ][@rank=%d] finalized MPI\n", rank);
//    return EXIT_SUCCESS;
////    return finalize_mpi_with_check(rank);
//  }
//  if (lm_index_max > num_of_wf_lm) {

  if ( lm_index_max >= num_of_wf_lm * (num_of_process + 1) ) {
//    std::cerr << "[ERROR] `lm_index_max` should not exceed `num_of_wf_lm`\n";
    fprintf(stderr, "[ERROR][@rank=%d] `lm_index_max` (=%d) should not exceed `num_of_wf_lm` (=%d)\n", rank, lm_index_max, num_of_wf_lm);
    return EXIT_FAILURE;
  } else {
    fprintf(stdout, "[ LOG ][@rank=%d] lm_index_start: %d / num_of_wf_to_read: %d\n", rank, lm_index_start, num_of_wf_to_read);
  }

  //// specify wf region to read
  long offset_lm;
  offset_lm = lm_index_start * bytes_per_wf;
  std::cout << "[@rank=" << rank << "]" << "[ LOG ] offset_lm = " << offset_lm << std::endl;


  //// Allocate memory for wavefunciton
  int num_of_elements_in_wf_stack = N_rho * num_of_wf_to_read;
  std::complex<double> *wf_read;
#ifdef HAVE_CUDA
  std::complex<double> *wf_read_aug;
  wf_read_aug = new std::complex<double>[num_of_elements_in_wf_stack + 2];
  wf_read_aug[0] = 0; wf_read_aug[num_of_elements_in_wf_stack + 1] = 0;  // set each ends as zeros
  wf_read = wf_read_aug + 1;
#else
//  if (num_of_wf_to_read != 0) {
  wf_read = new std::complex<double>[num_of_elements_in_wf_stack];
//  }
#endif

  //// Read wf data from file
#ifdef HAVE_MPI

  MPI_Datatype vec_type_wf_file;
  MPI_Type_vector(num_of_wf_to_read, N_rho, num_of_process, element_type, &vec_type_wf_file);
  MPI_Type_commit(&vec_type_wf_file);
  // Determine displacement (`disp`) for this process
//  MPI_Offset disp;
//MPI_Offset end_position;
//  disp = rank * num_of_wf_per_proc * element_type_size;
//  disp = num_of_wf_before_this_rank * element_type_size;
  // Declare variables for MPI writing process
  // Write tsurff quantities to file after setting view
  // - in this case, it corresponds to imposing a kind of mask to a file 
  //   so that a contiguous array of tsurff quantities can be stored at once 
  //   without looping over each time step
  // This is for `tsurffpsi.raw`
//  MPI_File_open(MPI_COMM_WORLD, "tsurffpsi.raw", MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &tsurff_psi_raw_file);
//  MPI_File_get_position(fh, &end_position);
  fprintf(stdout, "[ LOG ][@rank=%d] before MPI_File_set_view\n", rank);
  return_code = MPI_File_set_view(fh, offset_lm, element_type, vec_type_wf_file, "native", MPI_INFO_NULL);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_set_view"); }

  // The offset for `MPI_File_read_at()` is zero due to the file view
  fprintf(stdout, "[ LOG ][@rank=%d] before MPI_File_read\n", rank);
  MPI_Status read_status;
  return_code = MPI_File_read(fh, wf_read, num_of_elements_in_wf_stack, element_type, &read_status);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_read"); }

  fprintf(stdout, "[ LOG ][@rank=%d] after MPI_File_read\n", rank);

//  MPI_Status read_status;
//  return_code = MPI_File_read_at(fh, offset_lm, wf_read, num_of_wf_to_read * N_rho, MPI::DOUBLE_COMPLEX, &read_status);
//  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_read_at"); }
#else
//  fh.read((char *) wf_read, N_rho * num_of_wf_to_read * sizeof(std::complex<double>));
  fh.read((char *) wf_read, num_of_elements_in_wf_stack * sizeof(std::complex<double>));
  if (!fh) { 
    fprintf(stderr, "[ERROR] only `%ld` elements could be read\n", fh.gcount());
    fh.close();
    return EXIT_FAILURE;
  }
#endif


  //// Close file
#ifdef HAVE_MPI
  return_code = MPI_File_close(&fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_close"); }
#else
  fh.close();
#endif



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
  //
  // declare pointers for l-indenpendent, double arrays && allocation
  rho_array = new double[N_rho];
  scalarpot_array = new double[N_rho];
  imagpot_array = new double[N_rho];

  // fill out values
  double rho_value;
  for (i=0; i<N_rho; i++) {
    rho_value = delta_rho * (i+1);
    rho_array[i] = rho_value;
    scalarpot_array[i] = atomic_potential(rho_value, 0, 0, 0, 0);
    imagpot_array[i] = imaginarypot(i, 0, 0, 0, N_rho);
  }



  fprintf(stdout, "[ LOG ][@rank=%d] before evaluating tridiags\n", rank);
  
  //// Allocate memory for storing propagators for each `l` values
  std::complex<double> *tridiags_unitary_stack[NUM_OF_ARRAY_IN_TRIDIAGS], *tridiags_unitary_inv_stack[NUM_OF_ARRAY_IN_TRIDIAGS];
  for (i=0; i<NUM_OF_ARRAY_IN_TRIDIAGS; ++i) {

//    if (num_of_wf_to_read == 0) { break; }

    tridiags_unitary_stack[i] = new std::complex<double>[N_rho * num_of_wf_to_read];
    if (tridiags_unitary_stack[i] == NULL) { return error_and_exit(rank, EXIT_FAILURE, "malloc"); }
    tridiags_unitary_inv_stack[i] = new std::complex<double>[N_rho * num_of_wf_to_read];
    if (tridiags_unitary_inv_stack[i] == NULL) { return error_and_exit(rank, EXIT_FAILURE, "malloc"); }
  }

  fprintf(stdout, "[ LOG ][@rank=%d] right before evaluating tridiags\n", rank);

  //// Store propagators to each stacks
  // Prepare some variables
  long lm_index_offset;
  std::complex<double> *tridiags_unitary[NUM_OF_ARRAY_IN_TRIDIAGS], *tridiags_unitary_inv[NUM_OF_ARRAY_IN_TRIDIAGS];
  // Start storing
  int lm_index_index;
  for (lm_index=lm_index_start, lm_index_index=0; lm_index<lm_index_max; lm_index += num_of_process, ++lm_index_index) {

//    if (num_of_wf_to_read == 0) { break; }

    fprintf(stdout, "[ LOG ][@rank=%d] at the head lm_index loop, lm_index: %d\n", rank, lm_index);

    //// Retrieve `l` and `m` quantum numbers
    if (0 != get_ell_and_m_from_lm_index(lm_index, &l, &m, para_ini.getLong("initial-m"), g_prop.dimens())) { 
      fprintf(stderr, "[ERROR] during `get_ell_and_m_from_lm_index`\n");
      return EXIT_FAILURE; 
    }

    
    //// Evaulate unitary propagator
    lm_index_offset = lm_index - lm_index_start;

    // sign = -1 for explicit half time propagation
    for (i=0; i<NUM_OF_ARRAY_IN_TRIDIAGS; ++i) { 
      tridiags_unitary[i] = tridiags_unitary_stack[i] + lm_index_index * N_rho; 
//      tridiags_unitary[i] = tridiags_unitary_stack[i] + lm_index_offset * N_rho; 
    }
    sign = -1;
    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis_simple(
        l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,tridiags_unitary);

    // sign = 1 for implicit half time propagation
    for (i=0; i<NUM_OF_ARRAY_IN_TRIDIAGS; ++i) { 
      tridiags_unitary_inv[i] = tridiags_unitary_inv_stack[i] + lm_index_index * N_rho; 
//      tridiags_unitary_inv[i] = tridiags_unitary_inv_stack[i] + lm_index_offset * N_rho; 
    }
    sign = 1;
    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis_simple(
        l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,tridiags_unitary_inv);
    
    fprintf(stdout, "[ LOG ][@rank=%d] at the tail lm_index loop, lm_index: %d\n", rank, lm_index);

  } 


  fprintf(stdout, "[ LOG ][@rank=%d] all unitary propagator stored\n", rank);


  
  //// Set time index range for propagation
  long start_time_index = 0;
  double post_prop_duration = para_tsurff.getDouble("R-tsurff") / para_tsurff.getDouble("p-min-tsurff");
  double p_min_tsurff;
  long num_of_time_steps;
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
  long time_index_max = start_time_index + num_of_time_steps;
  long time_index;
  long num_of_steps_to_print_progress = 200; // [NOTE] to become global default config
  if (rank==0) {
    std::cout << "[@rank=" << rank << "]" << "[ LOG ] num_of_time_steps: " << num_of_time_steps << std::endl;
    std::cout << "[@rank=" << rank << "]" << "[ LOG ] post_prop_duration: " << post_prop_duration << std::endl;
  }

  //// Prepare memory for storing tsurff-quantity
  long tsurff_buffer_length = num_of_wf_to_read * num_of_time_steps;
  std::complex<double> *psi_R_arr = NULL;
  std::complex<double> *dpsi_drho_R_arr = NULL;
//  if (num_of_wf_to_read != 0) {
  psi_R_arr = new std::complex<double>[tsurff_buffer_length];
  dpsi_drho_R_arr = new std::complex<double>[tsurff_buffer_length];
//  }
  long index_at_R = g_prop.rindex(para_tsurff.getDouble("R-tsurff"));




  std::chrono::_V2::system_clock::time_point prop_start_time, prop_end_time;
  std::chrono::duration<double> elapsed_time_prop_total;
//  int elapsed_time_prop_total = -1;
//#ifdef HAVE_BOOST
//  boost::timer tim;
//#endif



  //// Propagate with tsurff quantity evaluation
  int return_code_prop = EXIT_FAILURE;
#ifdef HAVE_CUDA
  // forward tridiagonal multiplication configuration
  int num_of_thread_per_block = 128;
  int num_of_blocks_max = 32;
  int num_of_blocks = MIN((N_rho+num_of_thread_per_block-1)/num_of_thread_per_block, num_of_blocks_max);

  // Define grid and block dimension
  int block_dim3_in[3] = { num_of_thread_per_block, 1, 1 }, 
      grid_dim3_in[3] = { num_of_blocks, 1, 1 };

  int batch_stride = N_rho;


//#ifdef HAVE_BOOST
//    tim.restart();
//#endif 
  prop_start_time = std::chrono::high_resolution_clock::now();

  // running with gpu
  return_code_prop = cu_crank_nicolson_with_tsurff (
    index_at_R, delta_rho, start_time_index, num_of_time_steps,
    wf_read_aug, num_of_wf_to_read,
    psi_R_arr, dpsi_drho_R_arr, N_rho,
    tridiags_unitary_stack, tridiags_unitary_inv_stack, 
    num_of_steps_to_print_progress, rank,
    block_dim3_in, grid_dim3_in, batch_stride);


//#ifdef HAVE_BOOST
//  elapsed_time_prop_total = tim.elapsed();
//#endif // HAVE_BOOST
  prop_end_time = std::chrono::high_resolution_clock::now();

//  return_code_prop = tridiag_forward_backward (
//    N_rho, tridiags_unitary_stack, tridiags_unitary_inv_stack, wf_read_aug, 
//    start_time_index, time_index_max,
//    block_dim3_in, grid_dim3_in, 
//    num_of_wf_to_read, N_rho);
  if (return_code_prop != EXIT_SUCCESS) { 
    fprintf(stderr, "[ERROR] Abnormal exit from `cu_crank_nicolson_with_tsurff()`\n"); 
    return return_code_prop; 
  }
#else

//#ifdef HAVE_BOOST
//    tim.restart();
//#endif 
  prop_start_time = std::chrono::high_resolution_clock::now();

//  if (num_of_wf_to_read != 0) {
  //// running for non-gpu case
  return_code_prop = crank_nicolson_with_tsurff (
      index_at_R, delta_rho, start_time_index, num_of_time_steps, 
      wf_read, num_of_wf_to_read, psi_R_arr, dpsi_drho_R_arr, N_rho, 
      tridiags_unitary_stack, tridiags_unitary_inv_stack, 
      num_of_steps_to_print_progress, rank );
//  }

//#ifdef HAVE_BOOST
//  elapsed_time_prop_total = tim.elapsed();
//#endif // HAVE_BOOST
  prop_end_time = std::chrono::high_resolution_clock::now();

  if (return_code_prop != EXIT_SUCCESS) { return error_and_exit(rank, return_code_prop, "crank_nicolson_with_tsurff"); }
#endif // HAVE_CUDA


  //// Logging
  fprintf(stdout, "[ LOG ][@rank=%d] Propagation done\n", rank);

  elapsed_time_prop_total =  prop_end_time - prop_start_time;
  fprintf(stdout, "[ LOG ][@rank=%d] Elapsed time for propagation: `%f` seconds\n", rank, elapsed_time_prop_total.count());
//#ifdef HAVE_BOOST
//  cout << "time step took " << elapsed_time_prop_total.count() << " seconds" << endl;
//  cout << "time step took " << tim.elapsed() << " seconds" << endl;
//#endif


  //// Free arrays right after the propagation
  for (i=0; i<NUM_OF_ARRAY_IN_TRIDIAGS; ++i) {
    free(tridiags_unitary_stack[i]);
    free(tridiags_unitary_inv_stack[i]);
  }
  free(rho_array);
  free(scalarpot_array);
  free(imagpot_array);




  //// Write tsurff-quantities to files 
  string tsurff_psi_raw_file_name("tsurffpsi.raw"), tsurff_dpsidr_raw_file_name("tsurff-dpsidr.raw");
#ifdef HAVE_MPI

//  if (num_of_wf_to_read != 0) {

  // [DEPRECIATED] Define a MPI data type
  // : an array of blocks corresponding to each process, 
  // : thus separated with a certain interval
  // - The number of blocks is `num_of_time_steps`.
  // - The blocks are separated by a interval consisting of several `element_type`s.
  //   The length ot that interval is `num_of_wf_lm`.
  // - Each block consists of several `element_type`s.
  //   The number of the `element_type`s per block is `num_of_wf_to_read`
  
  
  // Define vector type for describing data strucutre for each time step
  MPI_Datatype vec_type_per_time_step;
  MPI_Type_vector(num_of_wf_to_read, 1, num_of_process, element_type, &vec_type_per_time_step);
  MPI_Type_commit(&vec_type_per_time_step);
  
  
  MPI_Datatype block_type;
  MPI_Type_contiguous(num_of_time_steps, vec_type_per_time_step, &block_type);
//  MPI_Type_vector(num_of_time_steps, num_of_wf_to_read, num_of_wf_lm, element_type, &block_type);
  MPI_Type_commit(&block_type);
  // Get `element_type`s size
//  int element_type_size;
//  MPI_Type_size(element_type, &element_type_size);
  // Determine displacement (`disp`) for this process
  MPI_Offset disp;
  MPI_Offset end_position;
//  disp = rank * num_of_wf_per_proc * element_type_size;
//  disp = num_of_wf_before_this_rank * element_type_size;
  disp = rank * element_type_size;

  // Declare variables for MPI writing process
  MPI_Status write_status;
  MPI_File tsurff_psi_raw_file, tsurff_dpsidr_raw_file;
  // Write tsurff quantities to file after setting view
  // - in this case, it corresponds to imposing a kind of mask to a file 
  //   so that a contiguous array of tsurff quantities can be stored at once 
  //   without looping over each time step
  // This is for `tsurffpsi.raw`
  MPI_File_open(working_world, "tsurffpsi.raw", MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &tsurff_psi_raw_file);
//  MPI_File_open(MPI_COMM_WORLD, "tsurffpsi.raw", MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &tsurff_psi_raw_file);
  MPI_File_get_position(tsurff_psi_raw_file, &end_position);
  MPI_File_set_view(tsurff_psi_raw_file, end_position+disp, element_type, block_type, "native", MPI_INFO_NULL);
  return_code = MPI_File_write(tsurff_psi_raw_file, psi_R_arr, tsurff_buffer_length, element_type, &write_status);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_write"); }
  return_code = MPI_File_close(&tsurff_psi_raw_file);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_close"); }
  // This is for `tsurff-dpsidr.raw`
  MPI_File_open(working_world, "tsurff-dpsidr.raw", MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &tsurff_dpsidr_raw_file);
//  MPI_File_open(MPI_COMM_WORLD, "tsurff-dpsidr.raw", MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &tsurff_dpsidr_raw_file);
  MPI_File_get_position(tsurff_dpsidr_raw_file, &end_position);
  MPI_File_set_view(tsurff_dpsidr_raw_file, end_position+disp, element_type, block_type, "native", MPI_INFO_NULL);
  return_code = MPI_File_write(tsurff_dpsidr_raw_file, dpsi_drho_R_arr, tsurff_buffer_length, element_type, &write_status);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_write"); }
  return_code = MPI_File_close(&tsurff_dpsidr_raw_file);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_close"); }

//  }

#else

  std::ofstream tsurff_psi_raw_file, tsurff_dpsidr_raw_file;
  // This is for `tsurffpsi.raw`
  tsurff_psi_raw_file.open(tsurff_psi_raw_file_name, std::ios::app | std::ios::binary);
  if (!tsurff_psi_raw_file.is_open()) { fprintf(stderr,"[ERROR] during opening file `%s`\n", tsurff_psi_raw_file_name.c_str()); }
  tsurff_psi_raw_file.write((char *) psi_R_arr, tsurff_buffer_length * sizeof(std::complex<double>));
  if(!tsurff_psi_raw_file || tsurff_psi_raw_file.fail()) { fprintf(stderr, "[ERROR] during writing `psi_R_arr` to file `tsurffpsi.raw`\n"); return EXIT_FAILURE; }
  else { fprintf(stdout, "[ LOG ] `%ld` elements has been written to file `%s`\n", tsurff_buffer_length, tsurff_psi_raw_file_name.c_str()); }
  tsurff_psi_raw_file.close();
  // This is for `tsurff-dpsidr.raw`
  tsurff_dpsidr_raw_file.open(tsurff_dpsidr_raw_file_name, std::ios::app | std::ios::binary);
  if (!tsurff_dpsidr_raw_file.is_open()) { fprintf(stderr,"[ERROR] during opening file `%s`\n", tsurff_dpsidr_raw_file_name.c_str()); }
  tsurff_dpsidr_raw_file.write((char *) dpsi_drho_R_arr, tsurff_buffer_length * sizeof(std::complex<double>));
  if(!tsurff_dpsidr_raw_file || tsurff_dpsidr_raw_file.fail()) { fprintf(stderr, "[ERROR] during writing `dpsi_drho_R_arr` to file `tsurff-dpsidr.raw`\n"); return EXIT_FAILURE; }
  else { fprintf(stdout, "[ LOG ] `%ld` elements has been written to file `%s`\n", tsurff_buffer_length, tsurff_dpsidr_raw_file_name.c_str()); }
  tsurff_dpsidr_raw_file.close();

#endif


  //// Write wf data to a file
#ifdef HAVE_MPI
//  if (num_of_wf_to_read != 0) {

  // Open file
  return_code = MPI_File_open(working_world, current_wf_bin_file_name.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
//  return_code = MPI_File_open(MPI_COMM_WORLD, current_wf_bin_file_name.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_open"); }
  // Write
  MPI_Status write_status_wf;
  return_code = MPI_File_write_at(fh, offset_lm, wf_read, num_of_wf_to_read * N_rho, element_type, &write_status_wf);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_write_at"); }

  // [NOTE] Check `write_status` to confirm how many elements has been written

  // Close file
  return_code = MPI_File_close(&fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_close"); }

//  }

#else
  fh.open(current_wf_bin_file_name.c_str(), std::ios::out | std::ios::binary);
  if (!fh.is_open()) { fprintf(stderr,"[ERROR] during opening file `%s`\n", current_wf_bin_file_name.c_str()); }
  fh.write((char *) wf_read, N_rho * num_of_wf_to_read * sizeof(std::complex<double>));
  if (!fh) { 
    fprintf(stderr, "[ERROR] during writing statefunction into file `%s`\n", current_wf_bin_file_name.c_str());
    fh.close();
    return EXIT_FAILURE;
  }
  fh.close();
#endif



  //// Free arrays
#ifdef HAVE_CUDA
  free(wf_read_aug);
#else
  free(wf_read);
#endif
  free(psi_R_arr);
  free(dpsi_drho_R_arr);





  //// MPI Finalization
#ifdef HAVE_MPI
  
//  MPI_Barrier(MPI_COMM_WORLD);

  
//fin:
  return finalize_mpi_with_check(rank);
//  return_code = MPI_Finalize();
//  if (return_code == MPI_SUCCESS) {
//    fprintf(stdout, "[@rank=%d][ LOG ] MPI program has been finalized.\n", rank);
//  } else {
//    fprintf(stderr, "[ERROR] Something got wrong during finalization.\n");
//    return EXIT_FAILURE;
//  }
#endif

  //// End this program
  return EXIT_SUCCESS;
}
