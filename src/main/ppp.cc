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

//// Define macro function
#define MIN(x,y) ((x<y)?x:y)


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
  int lm_index; // iteration for lm_unifired index
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



  //// Determine the number of basis wf (i.e. `wf_lm`) in the wf file
  long num_of_wf_lm = g_prop.num_of_phi_lm();
  if (num_of_wf_lm < 0) { std::cerr << "[ERROR] during `g_prop.num_of_phi_lm()`\n"; return EXIT_FAILURE; }
  if (rank == 0) { fprintf(stdout, "[ LOG ][@rank=%d] num_of_wf_lm = %ld\n", rank, num_of_wf_lm); }
  const int num_of_remaining_wf = num_of_wf_lm % num_of_process;



#ifdef HAVE_MPI

  //// Drop process that may be unused
  
  // Construct `world_group` representing the total set of processes
  MPI_Group world_group = MPI_GROUP_EMPTY;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  // Construct `working_group` representing a subset of processes which will actually run
  // .. In case where the number of process is larger than the num_of_wf_lm,
  // .. processes with last few ranks is dropped out and finalize before entering the calculation part
  MPI_Group working_group  = MPI_GROUP_EMPTY;
  int working_rank_range_incl[1][3] = { { 0, num_of_process-1, 1 } };
  if ( num_of_wf_lm < num_of_process ) { working_rank_range_incl[0][1] = num_of_remaining_wf-1; }
  MPI_Group_range_incl(world_group, 1, working_rank_range_incl, &working_group);

  // Construct MPI world `working_world` which consists of a group, the `working_group`
  // .. from now on, the `working_world` should be used instead of `MPI_COMM_WORLD`
  MPI_Comm working_world;
  MPI_Comm_create(MPI_COMM_WORLD, working_group, &working_world);

  // Finalize and exit the program if myself, represented by my `rank` does not need to run
  // .. most probably due to the lack of the number of radials wavefunctions to be calculated.
  if (working_world == MPI_COMM_NULL) { return finalize_mpi_with_check(rank); }

#endif // HAVE_MPI




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




  //// Determine the number of wf_lm to read and process: `num_of-wf_to_read` and the range of their `lm_index`s
  const int num_of_wf_per_proc_max = ( num_of_wf_lm + (num_of_process - 1) ) / num_of_process;
  const int num_of_wf_per_proc_min = num_of_wf_lm / num_of_process;
  const int num_of_wf_to_read = num_of_wf_per_proc_min + (rank < num_of_remaining_wf);
  const int lm_index_start = rank;
  const int lm_index_max = lm_index_start + num_of_wf_to_read * num_of_process;

  if ( (num_of_wf_lm < num_of_process) && (rank >= num_of_remaining_wf) ) {
    if ( (num_of_wf_to_read != 0) || (lm_index_max != lm_index_start) ) {
      fprintf(stderr, "[ERROR][@rank=%d] There is an inconsistency\n", rank);
      return EXIT_FAILURE;
    }
  }

  if ( lm_index_max >= num_of_wf_lm * (num_of_process + 1) ) {
    fprintf(stderr, "[ERROR][@rank=%d] `lm_index_max` (=%d) should not exceed `num_of_wf_lm` (=%ld)\n", rank, lm_index_max, num_of_wf_lm);
    return error_and_exit(rank, EXIT_FAILURE, "`lm_index_max` check");
  } else {
    fprintf(stdout, "[ LOG ][@rank=%d] lm_index_start: %d / num_of_wf_to_read: %d\n", rank, lm_index_start, num_of_wf_to_read);
  }




  //// Allocate memory for wavefunciton
  int num_of_elements_in_wf_stack = N_rho * num_of_wf_to_read;
  std::complex<double> *wf_read;
#ifdef HAVE_CUDA
  std::complex<double> *wf_read_aug;
  wf_read_aug = new std::complex<double>[num_of_elements_in_wf_stack + 2];
  wf_read_aug[0] = 0; wf_read_aug[num_of_elements_in_wf_stack + 1] = 0;  // set each ends as zeros
  wf_read = wf_read_aug + 1;
#else
  wf_read = new std::complex<double>[num_of_elements_in_wf_stack];
#endif



  //// Read wf data from file
#ifdef HAVE_MPI

  // specify wf region to read
  long offset_lm;
  offset_lm = lm_index_start * bytes_per_wf;
  std::cout << "[@rank=" << rank << "]" << "[ LOG ] offset_lm = " << offset_lm << std::endl;

  // Define MPI user-defined type for describing the wavefunction file
  MPI_Datatype vec_type_wf_file;
  MPI_Type_vector(num_of_wf_to_read, N_rho, num_of_process * N_rho, element_type, &vec_type_wf_file);
  MPI_Type_commit(&vec_type_wf_file);
  
  // Apply the defined MPI Datatype to the wavefunction file
  return_code = MPI_File_set_view(fh, offset_lm, element_type, vec_type_wf_file, "native", MPI_INFO_NULL);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_set_view"); }

  // [NOTE] The offset for `MPI_File_read_at()` is zero due to the file view
  MPI_Status read_status;
  return_code = MPI_File_read(fh, wf_read, num_of_elements_in_wf_stack, element_type, &read_status);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_read"); }
  
  // Get the number of elements that has been read from `MPI_File_read()`
  int num_of_elements_read = -1;
  return_code = MPI_Get_count(&read_status, element_type, &num_of_elements_read);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_Get_count"); }
  
  // Cehck whether all target elements has been successfully read from the target file
  if (num_of_elements_read != num_of_elements_in_wf_stack) {
    fprintf(stderr, "[ERROR][@rank=%d] Number of elements to be read (=%d) is different from the actual number of elements that has been read (=%d) from file `%s`\n", rank, num_of_elements_in_wf_stack, num_of_elements_read, current_wf_bin_file_name.c_str());
    return error_and_exit(rank, EXIT_FAILURE, "MPI_File_read");
  }


//  MPI_Status read_status;
//  return_code = MPI_File_read_at(fh, offset_lm, wf_read, num_of_wf_to_read * N_rho, MPI::DOUBLE_COMPLEX, &read_status);
//  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_read_at"); }
#else
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



  
  //// Allocate memory for storing propagators for each `l` values
  std::complex<double> *tridiags_unitary_stack[NUM_OF_ARRAY_IN_TRIDIAGS], *tridiags_unitary_inv_stack[NUM_OF_ARRAY_IN_TRIDIAGS];
  for (i=0; i<NUM_OF_ARRAY_IN_TRIDIAGS; ++i) {
    tridiags_unitary_stack[i] = new std::complex<double>[N_rho * num_of_wf_to_read];
    if (tridiags_unitary_stack[i] == NULL) { return error_and_exit(rank, EXIT_FAILURE, "malloc"); }
    tridiags_unitary_inv_stack[i] = new std::complex<double>[N_rho * num_of_wf_to_read];
    if (tridiags_unitary_inv_stack[i] == NULL) { return error_and_exit(rank, EXIT_FAILURE, "malloc"); }
  }


  //// Store propagators to each stacks
  // Prepare some variables
  long lm_index_offset;
  std::complex<double> *tridiags_unitary[NUM_OF_ARRAY_IN_TRIDIAGS], *tridiags_unitary_inv[NUM_OF_ARRAY_IN_TRIDIAGS];
  // Start storing
  int lm_index_index;
  for (lm_index=lm_index_start, lm_index_index=0; lm_index<lm_index_max; lm_index += num_of_process, ++lm_index_index) {
    
    //// Retrieve `l` and `m` quantum numbers
    if (0 != get_ell_and_m_from_lm_index(lm_index, &l, &m, para_ini.getLong("initial-m"), g_prop.dimens())) { 
      fprintf(stderr, "[ERROR][@rank=%d] during `get_ell_and_m_from_lm_index`\n", rank);
      return EXIT_FAILURE; 
    }
    
    //// Evaulate unitary propagator
    lm_index_offset = lm_index - lm_index_start;

    // sign = -1 for explicit half time propagation
    for (i=0; i<NUM_OF_ARRAY_IN_TRIDIAGS; ++i) { 
      tridiags_unitary[i] = tridiags_unitary_stack[i] + lm_index_index * N_rho; 
    }
    sign = -1;
    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis_simple(
        l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,tridiags_unitary);

    // sign = 1 for implicit half time propagation
    for (i=0; i<NUM_OF_ARRAY_IN_TRIDIAGS; ++i) { 
      tridiags_unitary_inv[i] = tridiags_unitary_inv_stack[i] + lm_index_index * N_rho; 
    }
    sign = 1;
    evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis_simple(
        l,sign,N_rho,Z,delta_rho,delta_t,rho_array,scalarpot_array,imagpot_array,tridiags_unitary_inv);
    
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
      fprintf(stderr, "[ERROR][@rank=%d] either `num-of-time-steps` or `p-min-tsurff` should be set\n", rank);
      return EXIT_FAILURE;
    }
    post_prop_duration = para_tsurff.getDouble("R-tsurff") / p_min_tsurff;
    num_of_time_steps = long(post_prop_duration / delta_t); 
  }
  long time_index_max = start_time_index + num_of_time_steps;
  long time_index;
  long num_of_steps_to_print_progress = 200; // [NOTE] to become global default config
  if (rank==0) {
    fprintf(stdout, "[ LOG ][@rank=%d] num_of_time_steps: %ld\n", rank, num_of_time_steps);
    fprintf(stdout, "[ LOG ][@rank=%d] post_prop_duration: %.3f [a.u.]\n", rank, post_prop_duration);
  }



  //// Prepare memory for storing tsurff-quantity
  long tsurff_buffer_length = num_of_wf_to_read * num_of_time_steps;
  std::complex<double> *psi_R_arr = NULL;
  std::complex<double> *dpsi_drho_R_arr = NULL;
  psi_R_arr = new std::complex<double>[tsurff_buffer_length];
  dpsi_drho_R_arr = new std::complex<double>[tsurff_buffer_length];
  long index_at_R = g_prop.rindex(para_tsurff.getDouble("R-tsurff"));



  //// Prepare timer
  std::chrono::_V2::system_clock::time_point prop_start_time, prop_end_time;
  std::chrono::duration<double> elapsed_time_prop_total;



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


  prop_start_time = std::chrono::high_resolution_clock::now();

  return_code_prop = cu_crank_nicolson_with_tsurff (
    index_at_R, delta_rho, start_time_index, num_of_time_steps,
    wf_read_aug, num_of_wf_to_read,
    psi_R_arr, dpsi_drho_R_arr, N_rho,
    tridiags_unitary_stack, tridiags_unitary_inv_stack, 
    num_of_steps_to_print_progress, rank,
    block_dim3_in, grid_dim3_in, batch_stride);


  prop_end_time = std::chrono::high_resolution_clock::now();

  if (return_code_prop != EXIT_SUCCESS) { 
    fprintf(stderr, "[ERROR] Abnormal exit from `cu_crank_nicolson_with_tsurff()` with code: `%d`\n", return_code_prop);
    return EXIT_FAILURE;
  }

#else // without CUDA

  prop_start_time = std::chrono::high_resolution_clock::now();

  //// running for non-gpu case
  return_code_prop = crank_nicolson_with_tsurff (
      index_at_R, delta_rho, start_time_index, num_of_time_steps, 
      wf_read, num_of_wf_to_read, psi_R_arr, dpsi_drho_R_arr, N_rho, 
      tridiags_unitary_stack, tridiags_unitary_inv_stack, 
      num_of_steps_to_print_progress, rank );

  prop_end_time = std::chrono::high_resolution_clock::now();

  if (return_code_prop != EXIT_SUCCESS) { return error_and_exit(rank, return_code_prop, "crank_nicolson_with_tsurff"); }

#endif // HAVE_CUDA


  //// Logging
  fprintf(stdout, "[ LOG ][@rank=%d] Propagation done\n", rank);

  elapsed_time_prop_total =  prop_end_time - prop_start_time;
  fprintf(stdout, "[ LOG ][@rank=%d] Elapsed time for propagation: `%f` seconds\n", rank, elapsed_time_prop_total.count());



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
//  MPI_Type_vector(num_of_wf_per_proc_max, 1, num_of_process * 1, element_type, &vec_type_per_time_step);
  MPI_Type_vector(num_of_wf_to_read, 1, num_of_process * 1, element_type, &vec_type_per_time_step);
  MPI_Type_commit(&vec_type_per_time_step);
  
  
  MPI_Datatype block_type;
  MPI_Type_create_hvector(num_of_time_steps, 1, num_of_wf_lm * element_type_size, vec_type_per_time_step, &block_type);
//  MPI_Type_contiguous(num_of_time_steps, vec_type_per_time_step, &block_type);
//  MPI_Type_vector(num_of_time_steps, num_of_wf_to_read, num_of_wf_lm, element_type, &block_type);
  MPI_Type_commit(&block_type);
  // Determine displacement (`disp`) for this process
  MPI_Offset disp;
  MPI_Offset end_position;
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


#else // without MPI

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
  


  // Define MPI user-defined type for describing the wavefunction file
//  MPI_Datatype vec_type_wf_file;
//  MPI_Type_vector(num_of_wf_to_read, N_rho, num_of_process * N_rho, element_type, &vec_type_wf_file);
//  MPI_Type_commit(&vec_type_wf_file);
  



  // Open file
  return_code = MPI_File_open(working_world, current_wf_bin_file_name.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
//  return_code = MPI_File_open(MPI_COMM_WORLD, current_wf_bin_file_name.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_open"); }
  
  // Apply the defined MPI Datatype to the wavefunction file
  return_code = MPI_File_set_view(fh, offset_lm, element_type, vec_type_wf_file, "native", MPI_INFO_NULL);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_set_view"); }

  // Write
  MPI_Status write_status_wf;
  // [NOTE] The offset for `MPI_File_write_at()` is zero due to the file view
  return_code = MPI_File_write_at(fh, 0, wf_read, num_of_wf_to_read * N_rho, element_type, &write_status_wf);
//  return_code = MPI_File_write_at(fh, offset_lm, wf_read, num_of_wf_to_read * N_rho, element_type, &write_status_wf);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_write_at"); }

  int num_of_elements_written = -1;
  return_code = MPI_Get_count(&write_status_wf, element_type, &num_of_elements_written);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_Get_count"); }

  if (num_of_elements_written != num_of_wf_to_read * N_rho) {
    fprintf(stderr, "[ERROR][@rank=%d] Inconsistent writing file `%s`\n", rank, current_wf_bin_file_name.c_str());
    return error_and_exit(rank, EXIT_FAILURE, "MPI_File_Write_at");
  }

  // [NOTE] Check `write_status` to confirm how many elements has been written

  // Close file
  return_code = MPI_File_close(&fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_close"); }

#else // wihtout MPI

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



  //// Free memory allocated for arrays
#ifdef HAVE_CUDA
  free(wf_read_aug);
#else
  free(wf_read);
#endif
  free(psi_R_arr);
  free(dpsi_drho_R_arr);



  //// MPI Finalization
#ifdef HAVE_MPI
  return finalize_mpi_with_check(rank);
#endif



  //// End this program
  return EXIT_SUCCESS;

}
