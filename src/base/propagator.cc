#include "propagator.hh"

void evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis_simple(long l, long sign, long N_rho, double Z, double delta_rho, double delta_t, double *rho_array, double *scalarpot, double *imagpot, std::complex<double> *unitary_tridiags[]) {
  
  //// Get handles from tridiags, which is an array of 1D array pointers
  std::complex<double> *diag_unitary = unitary_tridiags[i_d];
  std::complex<double> *lower_offdiag_unitary = unitary_tridiags[i_ld] + 1;
  std::complex<double> *upper_offdiag_unitary = unitary_tridiags[i_ud];

  //// Set offdiag's boundary conditions
  lower_offdiag_unitary[-1] = 0.0;
  upper_offdiag_unitary[N_rho-1] = 0.0;

  //// Evaluate propagator
  evaluate_Numerov_boosted_CN_propagator_tridiags_for_sph_harm_basis(
      l, sign, N_rho, Z, delta_rho, delta_t, rho_array, scalarpot, imagpot, 
      diag_unitary, lower_offdiag_unitary, upper_offdiag_unitary );
}



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



int crank_nicolson_with_tsurff(
    int index_at_R, double delta_rho, int start_time_index, int num_of_time_steps, 
    std::complex<double> *wf_read, int num_of_wf_lm, 
    std::complex<double> *psi_R_arr, std::complex<double> *dpsi_drho_R_arr, long N_rho, 
    std::complex<double> *tridiags_unitary_stack[], std::complex<double> *tridiags_unitary_inv_stack[], 
    int num_of_steps_to_print_progress, int rank) 
{
  

  //// Declare/Define intermediate variables
  
  // time varialbes
  int time_index;
  int time_index_max = start_time_index + num_of_time_steps;
  long num_of_steps_done_so_far, time_index_from_zero;

  // lm-related variables
  int lm_index_offset, lm_index;
  int num_of_element_in_wf_lm = N_rho;
  int total_num_of_element_in_given_wf = num_of_element_in_wf_lm * num_of_wf_lm;
  long num_of_numbers_before_this_wf_lm;
  std::complex<double> *wf_lm = NULL, *wf_lm_mid = NULL;

  // tsurff quantity evaluation related variables
  const double one_over_12delta_rho = 1.0/(12.0*delta_rho);
  const double two_over_3delta_rho = 2.0/(3.0*delta_rho);
  int tsurff_buf_index;
  


  //// Start iteration over each `wf_lm` for `lm_index`
  
  // allocate memory for an intermediate array
  wf_lm_mid = new std::complex<double>[total_num_of_element_in_given_wf];

  // iterate over time
  for (time_index=start_time_index; time_index<time_index_max; ++time_index) {

    time_index_from_zero = time_index - start_time_index;

    // iterate over `wf_lm`
    for (lm_index_offset=0; lm_index_offset<num_of_wf_lm; ++lm_index_offset) {
  
      // Set wf_lm pointer
      num_of_numbers_before_this_wf_lm = lm_index_offset * N_rho;
      wf_lm = wf_read + num_of_numbers_before_this_wf_lm; 
  
      // Evaulate tsurff-related values
      tsurff_buf_index = time_index_from_zero * num_of_wf_lm + lm_index_offset;
      psi_R_arr[tsurff_buf_index] = wf_lm[index_at_R];
      dpsi_drho_R_arr[tsurff_buf_index] = two_over_3delta_rho*(wf_lm[index_at_R+1] - wf_lm[index_at_R-1]) - one_over_12delta_rho*(wf_lm[index_at_R+2]-wf_lm[index_at_R-2]);
  
    }


    //// Propagate
    // explicit half time propagation
    tridiag_mul_forward(tridiags_unitary_stack[i_ld], tridiags_unitary_stack[i_d], tridiags_unitary_stack[i_ud], wf_read, wf_lm_mid, total_num_of_element_in_given_wf);
    // implicit half time propagation
    tridiag_mul_backward(tridiags_unitary_inv_stack[i_ld], tridiags_unitary_inv_stack[i_d], tridiags_unitary_inv_stack[i_ud], wf_read, wf_lm_mid, total_num_of_element_in_given_wf);
  

    //// Logging
    if (((time_index + 1) % num_of_steps_to_print_progress) == 0) {
      num_of_steps_done_so_far = time_index_from_zero + 1;
      if (rank == 0) {
        fprintf(stdout, "[@rank=%d][ LOG ] num_of_steps_done_so_far = %ld / %d\n", 
            rank, num_of_steps_done_so_far, num_of_time_steps);
      }
    }

  }
  
  //// Free memory from inside
  free(wf_lm_mid);
  
  //// END OF ROUTINE
  return EXIT_SUCCESS;
}




