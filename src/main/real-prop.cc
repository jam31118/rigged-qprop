#include "real-prop.h"

#ifdef BOHM
#include "log.hh"
#include "lm.hh"
#include "velocity.hh"
#endif // BOHM

using std::endl;
using std::cout;

void print_banner() {
  fprintf(stdout, " IONIZATION\n");
  fprintf(stdout, " (C) Copyright by Bauer D and Koval P, Heidelberg (2005)\n");
  fprintf(stdout, " -------------------------------------------------------\n");
};


int real_prop(int argc, char **argv) {

  grid g_prop;
  wavefunction staticpot, wf, wf_initial;

  print_banner();

  // Set verbosity of stdout
  int iv            = 1;


  // Instantiate parameter objects from which several user-defined parameters 
  // will be extracted and used in this prgoram
  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");

  parameterListe default_param = get_default_parameter_list_object(); 


  // Get some parameters from parameter files
  const int my_m_quantum_num = para_ini.getLong("initial-m");
  const double nuclear_charge = para_ini.getDouble("nuclear-charge");

  long start_time_index;
  try { start_time_index = para_prop.getLong("start-time-index"); }
  catch (std::exception&) { start_time_index = 0; }
  

  // Output configuration
  // The log files are written with this `default_mode`
  // it is either "w" (new write; erasing previous file contents if any) or "a" (append mode)
  string default_mode_string;
  if (start_time_index == 0) { default_mode_string = string("w"); }
  else if (start_time_index > 0) { default_mode_string = string("a"); }
  else { std::cerr << "Unexpected start_time_index: " 
    << start_time_index << endl; return 1; }
  const char *default_mode = default_mode_string.c_str();


  // Declare and configure a grid for propagation
  const double delta_r = para_ini.getDouble("delta-r");
  g_prop.set_dim(para_ini.getLong("qprop-dim"));


  // Add a parameter which determines radial distance between R-tsurff and imag-pot
  double beyond_R_distance_temp;
  try { beyond_R_distance_temp = para_tsurff.getDouble("beyond-R"); }
  catch (std::exception&) { beyond_R_distance_temp = 0.0; }
  const double beyond_R_distance = beyond_R_distance_temp;
  double grid_size = para_tsurff.getDouble("R-tsurff") + beyond_R_distance + para_prop.getDouble("imag-width");
  
  g_prop.set_ngps(long(grid_size/delta_r), para_prop.getLong("ell-grid-size"), 1); 
  g_prop.set_delt(delta_r);
  g_prop.set_offs(0, 0, 0);
  
  
  // Construct Vector potential objects
  superposed_vecpot vecpot_x, vecpot_y, vecpot_z;
  construct_vecpot(g_prop.dimens(), para_prop, vecpot_x, vecpot_y, vecpot_z);


  // Confirm whether the vecpot object had been generated well.
  if ( g_prop.dimens() == 34 ) {
    print_superposed_vecpot(vecpot_z, "superposed_vecpot_z");
  } else if ( g_prop.dimens() == 44 ) {
    print_superposed_vecpot(vecpot_x, "superposed_vecpot_x");
    print_superposed_vecpot(vecpot_y, "superposed_vecpot_y");
  }


  // Determine duration of the vector potential
  // [NOTE] This should be determined by qprop-dimension and combination of several vecpot
  double pulse_duration;
  if ( g_prop.dimens() == 34 ) { 
    if ( vecpot_z.get_start_time() != 0 ) { std::cerr << "[ERROR] global time should be zero\n"; }
    pulse_duration = vecpot_z.get_duration(); 
  }
  else if ( g_prop.dimens() == 44 ) {
    if ( ( vecpot_x.get_start_time() != 0 ) || ( vecpot_y.get_start_time() != 0 ) ) { 
      std::cerr << "[ERROR] global start time should should be zero\n"; }
    double pulse_duration_x = vecpot_x.get_duration();
    double pulse_duration_y = vecpot_y.get_duration();
    if ( pulse_duration_x > pulse_duration_y ) { pulse_duration = pulse_duration_x; }
    else { pulse_duration = pulse_duration_y; }
  }


  //// Get and set PPP configuation
  bool use_ppp = false;
  try { use_ppp = para_tsurff.getBool("use-ppp"); }
  catch (std::exception&) {}

  // How long do the slowest electrons have time to reach the t-SURFF boundary
  // [NOTE] Keep in mind that the electron doesn't have a position 
  // and there's no concept like 'reach somewhere' for the electron especially in this quantum treatment.
  const double time_surff=para_tsurff.getDouble("R-tsurff")/para_tsurff.getDouble("p-min-tsurff");
  double duration = pulse_duration;
  if (!use_ppp) { duration += time_surff; }
  cout << "[ LOG ] pulse_duration : " << pulse_duration << endl;
  cout << "[ LOG ] time_surff : " << time_surff << endl;
  if (use_ppp) {
    cout << "[ LOG ] Run `PPP` for `time_surff`\n";
  }
  
  
  // The output files that will be created by `real-prop`
  string common_prefix("real-prop");
  string str_fname_logfi=common_prefix+string(".log");
  FILE* file_logfi = fopen_with_check(str_fname_logfi, default_mode);
  string str_fname_yield=common_prefix+string("-yield.dat");
  FILE* file_yield = fopen_with_check(str_fname_yield, default_mode);
  string str_fname_obser=common_prefix+string("-obser.dat");
  FILE* file_obser = fopen_with_check(str_fname_obser, default_mode);
  
  if (iv!=0) {
    fprintf(stdout, "%s will be (re)written.\n", str_fname_yield.c_str());
    fprintf(stdout, "%s will be (re)written.\n", str_fname_logfi.c_str());
  };


  // create an instance of the class for doing the tsurff related work
  tsurffSaveWF tsurff_save_wf(para_ini, para_prop, para_tsurff, g_prop);

  // the absorbing imaginary potential
  const long imag_potential_width=long(para_prop.getDouble("imag-width")/delta_r);
  imagpot imaginarypot(imag_potential_width);
  if ( print_imagpot(-1, NULL) != 0) { fprintf(stderr, "[ERROR] Failed to print imagpot\n"); };

  // set the binding potential and the hamiltonian
//  double alpha;
//  try { alpha = para_ini.getDouble("effpot-alpha"); }
//  catch (std::exception&) { alpha = 0.0; } // default value - if alpha = 0.0 means no effective potential but just become coulumb potential.
  scalarpot scalarpotx(para_ini.getDouble("nuclear-charge"), para_ini.getDouble("pot-cutoff"), get_effpot_alpha(para_ini));
  if ( print_scalarpot(-1, NULL) != 0) { fprintf(stderr, "[ERROR] Failed to print scalarpot\n"); };
  

  // Define hamiltonian object
  hamop hamilton;
  hamilton.init(g_prop, vecpot_x, vecpot_y, vecpot_z, scalarpotx, always_zero5, always_zero5, imaginarypot, always_zero2);


  // this is the linear and constant part of the Hamiltonian
  staticpot.init(g_prop.size()); 
  staticpot.calculate_staticpot(g_prop, hamilton);


  // Initialize wavefunction array 
  wf.init(g_prop.size()); 
  wf_initial.init(g_prop.size());


  // Load initial wavefunction from file
  string initial_wf_file_name = default_param.getString("initial_wf_file_name"); // string("ini-wf.bin"); // [NOTE] global variable -> 180-716 done.
  std::ifstream initial_wf_file(initial_wf_file_name, std::ios::binary);
  if ( !initial_wf_file.is_open() ) {
    std::cerr << "[ERROR] Failed to read file: " << initial_wf_file_name << endl;
    return 1; }
  if ( wf_initial.load_from_binary(initial_wf_file) ) {
    std::cerr << "[ERROR] Failed to load initial wavefunction from file with name: " 
      << initial_wf_file_name << endl;
    return 1; }
  initial_wf_file.close();
  

  // Load start (may be different from 'initial') wavefunction from file with name: 'start_wf_file_name'
  // [NOTE] 'start wavefunction' is different from 'initial wavefunction' if 'start_time_index' is not zero,
  // i.e. if the calculation doesn't start from zero (initial) time.
  string start_wf_file_name;
  if (start_time_index == 0) { start_wf_file_name = initial_wf_file_name; }
//  else if (start_time_index > 0) { start_wf_file_name = string("current-wf.bin"); } // [NOTE] global variable
  else if (start_time_index > 0) { start_wf_file_name = default_param.getString("current_wf_bin_file_name"); } // [NOTE] global variable -> 180716 done.
  else {
    fprintf(stderr, "[ERROR] Unexpected 'start_time_index': %ld\n", start_time_index);
    return 1; }
  std::ifstream start_wf_file(start_wf_file_name, std::ios::binary);
  if (!start_wf_file.is_open()) { 
    std::cerr << "[ERROR] Failed to read file: " << start_wf_file_name << endl; 
    return 1; }
  if (wf.load_from_binary(start_wf_file)) {
    std::cerr << "[ERROR] Failed to load starting wavefunction from file with name: " 
      << start_wf_file_name << endl;
    return 1; }
  start_wf_file.close();


  string str_fname_wf=common_prefix+string("-wf.dat");
  FILE* file_wf = fopen_with_check(str_fname_wf, "w");
  wf.dump_to_file_sh(g_prop, file_wf, 1); // initial wf is saved

  fprintf(stdout, "Norm of orbital: %le\n", wf.norm(g_prop));

  const double real_timestep = para_prop.getDouble("delta-t");
  long lno_of_ts = long( duration*1.0/real_timestep ) + 1;
  //
  // write to log file
  //
  fprintf(file_logfi, "Real-time propagation\n");
  fprintf(file_logfi, "Grid: \n");
  fprintf(file_logfi, "g_prop.ngps_x() = %ld\n", g_prop.ngps_x());
  fprintf(file_logfi, "g_prop.ngps_y() = %ld\n", g_prop.ngps_y());
  fprintf(file_logfi, "g_prop.ngps_z() = %ld\n", g_prop.ngps_z());
  fprintf(file_logfi, "g_prop.dimens() = %d\n\n", g_prop.dimens());
  fprintf(file_logfi, "g_prop.delt_x() = %15.10le\n", g_prop.delt_x());

  fprintf(file_logfi, "real_timestep     = %15.10le\n", real_timestep);
  fprintf(file_logfi, "lno_of_ts         = %ld\n", lno_of_ts);
  fprintf(file_logfi, "nuclear_charge    = %15.10le\n", nuclear_charge);
  fprintf(file_logfi, "str_fname_obser  = %s\n", str_fname_obser.c_str());
  fprintf(file_logfi, "str_fname_wf = %s\n", str_fname_wf.c_str());
  fprintf(file_logfi, "duration = %15.10le\n", duration);
  fflush(file_logfi);
  fclose(file_logfi);

  long ldumpwidth(1.0/real_timestep);  // output once every a.u. of time
  int me = 0; // dummy here


  // write vector potential to file
  string str_fname_vpot;
  if (start_time_index == 0) {
    str_fname_vpot = common_prefix+string("-vecpot.dat");
    FILE* file_vpot = fopen_with_check(str_fname_vpot, "w");
    double time = 0.0;
    if (g_prop.dimens() == 34) {
      for (long ts=0; ts<lno_of_ts; ts++) {
        if (ts%ldumpwidth==0)
          { fprintf(file_vpot, "%15.10le %15.10le\n", 
              time, vecpot_z(time, me)); }
        time += real_timestep;
      }
    } else if (g_prop.dimens() == 44) {
      for (long ts=0; ts<lno_of_ts; ts++) {
        if (ts%ldumpwidth==0)
          { fprintf(file_vpot, "%15.10le %15.10le %15.10le\n", 
              time, vecpot_x(time, me), vecpot_y(time, me)); }
        time += real_timestep;
      }
    } else { 
      std::cerr << "[ERROR] Unexpected grid dimension: " << g_prop.dimens() << endl;
      return 1; 
    }
    fclose(file_vpot);
  }


  // Configure saving wavefunction
  string current_wf_bin_file_name = string("current-wf.bin");


  // Configure Backup
  long wf_backup_period_temp;
  try { wf_backup_period_temp = para_ini.getLong("wf-backup-period"); }
  catch (std::exception&) { wf_backup_period_temp = 2000; } // default values
  if (wf_backup_period_temp < 0) { 
    std::cerr << "[ERROR] wavefunction backup period should be positive integer.\n"
      << "[ LOG ] Input backup period = " << wf_backup_period_temp << endl;
    return 1; }
  const long wf_backup_period = wf_backup_period_temp;
  

  // ********************************************************
  // ***** real time propagation ****************************
  // ********************************************************
  cplxd timestep=cplxd(real_timestep, 0.0);
  cplxd P;
  double N;
  
  // get 'num-time-index' from parameter file
  long num_time_index, max_time_index;
  try { num_time_index = para_prop.getLong("num-time-index"); }
  catch (std::exception&) {
    max_time_index = lno_of_ts; 
    num_time_index = max_time_index - start_time_index;
  }
  max_time_index = start_time_index + num_time_index;
  if (max_time_index > lno_of_ts) {
    std::cerr << "[ERROR] max_time_index " << max_time_index
      << "exceeds total number of timesteps " << lno_of_ts << endl;
    std::cerr << "[ LOG ] Aborting . . ." << endl;
    return 1;
  }

#ifdef HAVE_BOOST
  boost::timer tim;
#endif
#ifdef HAVE_BOOST
    tim.restart();
#endif




#ifdef BOHM
  
  //// Configuration
  const int N_s = 4, N_rho = g_prop.ngps_x();
  const int N_r_dim = 3, N_l = g_prop.ngps_y(), qprop_dim=g_prop.dimens();
  int N_p_;
  try { N_p_ = para_prop.getLong("num-of-particle"); }
  catch (std::exception&) { return return_with_mesg("Please specify the 'num-of-particle'"); }
  const int N_p = N_p_;
  if (N_p <= 0) { return return_with_mesg("The number of particle should be positive"); }

  // N_l + 1; // initialized with forbidden m value
  const int initial_m = my_m_quantum_num; 
  const double rho_max = 5.0;

  const int N_t = num_time_index + 1;
//  const double delta_t = 0.2;
//  const double t_0 = 0.2;

  int N_lm;
  if ( eval_N_lm(qprop_dim, N_l, &N_lm) != EXIT_SUCCESS) 
  { return return_with_mesg("Failed to evaluate 'N_lm'"); }
  if ( N_lm != g_prop.num_of_phi_lm() ) 
  { return return_with_mesg("Inconsistent N_lm"); }

  int *l_arr = new int[N_lm];
  int *m_arr = new int[N_lm];
  if ( eval_l_m_arr(qprop_dim, N_l, l_arr, m_arr, initial_m) != EXIT_SUCCESS )
  { return return_with_mesg("Failed to evaluate 'l_arr', 'm_arr'"); }

  
  
  
  //// Define useful variables
  int return_code = EXIT_FAILURE;


  //// Allocate memory for data

  // Some variable
  const int r_p_arr_size = N_r_dim * N_p;

  // `r_p_arr`
  double **r_p_arr = new double*[N_r_dim];
  double *r_p_arr_1d = new double[r_p_arr_size];
  for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
    r_p_arr[i_r_dim] = r_p_arr_1d + i_r_dim * N_p;
  }

  // `v_p_arr`
  double **v_p_arr = new double*[N_r_dim];
  double *v_p_arr_1d = new double[r_p_arr_size];
  for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
    v_p_arr[i_r_dim] = v_p_arr_1d + i_r_dim * N_p;
  }
  
  // `psi_arr`
  cplxd **psi_arr = new cplxd*[N_lm];
  psi_arr[0] = wf.begin();  // new cplxd[N_lm*N_rho];
  for (int i_lm=0; i_lm<N_lm; i_lm++) 
  { psi_arr[i_lm] = psi_arr[0] + i_lm * N_rho; }
  
  // `r_p_t_arr`
  double **r_p_t_arr = new double*[N_t];
  double *r_p_t_arr_1d = new double[N_t*r_p_arr_size];
  for (int i_t = 0; i_t < N_t; i_t++) 
  { r_p_t_arr[i_t] = r_p_t_arr_1d + i_t * r_p_arr_size; }
  
  // `v_p_t_arr`
  double **v_p_t_arr = new double*[N_t];
  double *v_p_t_arr_1d = new double[N_t*r_p_arr_size];
  for (int i_t = 0; i_t < N_t; i_t++) 
  { v_p_t_arr[i_t] = v_p_t_arr_1d + i_t * r_p_arr_size; }


  // Define variables for looping
  double *r_p_p, *r_p_p_max = r_p_p_max = r_p_arr_1d + r_p_arr_size;
  



  //// Initialization

  const double rho_p_lim[2] = { 0.0, para_tsurff.getDouble("R-tsurff") };

  // `rho_p_arr`
  double *rho_p_arr = new double[N_p];
  double delta_rho_p = (20.0 - rho_p_lim[0]) / (N_p+1);
  for (int i_p = 0; i_p < N_p; i_p++) {
    rho_p_arr[i_p] = delta_rho_p * (i_p + 1);
  }
//  double rho_p_arr[N_p] = {0.5, 1.0, 2.1, 2.4};
  for (int i_p = 0; i_p < N_p; i_p++) {
    r_p_arr[0][i_p] = rho_p_arr[i_p];
    r_p_arr[1][i_p] = 0.3 * M_PI;
    r_p_arr[2][i_p] = 0.0;
  }

  // Construct radial coordinate grid points `rho_arr`
//  const double delta_rho = (rho_max - 0) / (N_rho);
  const double delta_rho = delta_r;
 
  double *rho_arr = new double[N_rho];
  for (int i_rho = 0; i_rho < N_rho; i_rho++) 
  { rho_arr[i_rho] = g_prop.r(i_rho); } // (i_rho+1) * delta_rho; 

  // `v_p_arr`
  return_code = eval_v_p_arr_for_sph_harm_basis(
      N_s, N_p, N_rho, N_lm, r_p_arr, (const cplxd**) psi_arr, 
      rho_arr, l_arr, m_arr, rho_p_lim, v_p_arr);

  if (return_code != EXIT_SUCCESS) {
    fprintf(stderr, "[ERROR] Failed to run "
        "'eval_v_p_arr_for_sph_harm_basis()'");
    return return_code;
  }

  // `r_p_t_arr_1d` and `v_p_t_arr_1d` at t=t_0
  std::copy(r_p_arr_1d,r_p_arr_1d+r_p_arr_size,r_p_t_arr_1d+0*r_p_arr_size);
  std::copy(v_p_arr_1d,v_p_arr_1d+r_p_arr_size,v_p_t_arr_1d+0*r_p_arr_size);





  const int N_st = 4;

  double A[N_st][N_st] = { {0,0,0,0},{0.5,0,0,0},{0,0.5,0,0,},{0,0,1,0} };
  const double c[N_st] = {0.0, 0.5, 0.5, 1.0};
  const double b[N_st] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
  
  //// Allocate temporary buffer
  double ***k_st_arr = new double**[N_st];
  double **k_st_arr_2d = new double*[N_st*N_r_dim];
  double *k_st_arr_1d = new double[N_st*N_r_dim*N_p];
  for (int i_st = 0; i_st < N_st; i_st++) {
    k_st_arr[i_st] = k_st_arr_2d + i_st * N_r_dim;
    for (int i_st_r_dim = 0; i_st_r_dim < N_st*N_r_dim; i_st_r_dim++) {
      k_st_arr_2d[i_st_r_dim] = k_st_arr_1d + i_st_r_dim * N_p;
    }
  }



  double **r_p_i_st_arr = new double*[N_r_dim];
  double *r_p_i_st_arr_1d = new double[r_p_arr_size];
  for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
    r_p_i_st_arr[i_r_dim] = r_p_i_st_arr_1d + i_r_dim * N_p;
  }
  
  double *k_p; // *k_p_max;
  double a_i_k;



  double _delta_t, _delta_inner_t;
  int i_st, i_k, i_t = 1;

#endif // BOHM



  for (long ts=start_time_index; ts < max_time_index; ts++) {
    const double time=real_timestep*double(ts);

    if (ts%ldumpwidth==0) {
      // calculate total energy, projection onto initial state, norm, and <z>
      double E_tot = real(wf.energy(0.0, g_prop, hamilton, me, staticpot, nuclear_charge));
      // [NOTE] since the grid of both wavefunction are same, 
      // we can use this simple opeerator (`*`), instead of wavefunction::project() method
//      P = wf.project(g_prop, g_load, wf_load, 0);
      P = wf_initial * wf; 
      N = wf.norm(g_prop);
      double z_expect = real(wf.expect_z(g_prop));
      fprintf(file_obser, "%15.10le %15.10le %15.10le %15.10le %15.10le\n", 
          time, E_tot, real(conj(P)*P), N, z_expect);
    };
    
    // save the orbitals \varphi_{\ell}(\RI) and the derivative \partial_r\varphi_{\ell}(r)|_{r=\RI}
    tsurff_save_wf(wf);

#ifdef BOHM
    _delta_t = real_timestep;
    
    //// Eval k1 (assuming v_p_arr is already k1, provided c[0] == 0)
    i_st = 0;
    std::copy(v_p_arr_1d, v_p_arr_1d + r_p_arr_size, k_st_arr[i_st][0]);



    for (i_st = 1; i_st < N_st; i_st++) {
    
      _delta_inner_t = c[i_st-1] - c[i_st];

      // Eval y_i_st
      std::copy(r_p_arr_1d, r_p_arr_1d + r_p_arr_size, r_p_i_st_arr_1d);
      for (i_k=0; i_k < i_st; i_k++) 
      {
        a_i_k = A[i_st][i_k];

        for(r_p_p = r_p_i_st_arr_1d, 
            r_p_p_max = r_p_i_st_arr_1d + r_p_arr_size, 
            k_p = k_st_arr[i_k][0]; 
            r_p_p < r_p_p_max; r_p_p++, k_p++) 
        {
          *r_p_p += a_i_k * *k_p;
        }

      } // for-loop : i_k
      

      // Eval psi at t_i_st = t_n + c_i_st * _delta_t
      if (_delta_inner_t != 0) {
        wf.propagate(
            timestep, time, g_prop, hamilton, me, staticpot, 
            my_m_quantum_num, nuclear_charge);
//        return_code = propa_psi_arr(
//            psi_arr, _delta_inner_t*_delta_t, N_rho,N_lm,qprop_dim,initial_m);
//        if (return_code != EXIT_SUCCESS) {
//          return return_with_mesg("Failed during propagation of 'psi_arr'");
//        }
      }


      // Eval k
      return_code = eval_v_p_arr_for_sph_harm_basis(
          N_s, N_p, N_rho, N_lm, r_p_i_st_arr, (const cplxd **) psi_arr, 
          rho_arr, l_arr, m_arr, rho_p_lim, k_st_arr[i_st]);
      if (return_code != EXIT_SUCCESS) {
        return return_with_mesg("Failed to run 'eval_psi_and_dpsidx_arr()");
      }

    } // for-loop : i_st


    for (i_st = 0; i_st < N_st; i_st++) {
      for(r_p_p = r_p_arr_1d, 
          r_p_p_max = r_p_arr_1d + r_p_arr_size, 
          k_p = k_st_arr[i_st][0]; 
          r_p_p < r_p_p_max; r_p_p++, k_p++) 
      {
        *r_p_p += _delta_t * b[i_st] * *k_p;
      }
    } // for-loop : i_st

    //// Evaluate velocity vector for each particle
    return_code = eval_v_p_arr_for_sph_harm_basis(
      N_s, N_p, N_rho, N_lm, r_p_arr, (const cplxd**) psi_arr, 
        rho_arr, l_arr, m_arr, rho_p_lim, v_p_arr);
    if (return_code != EXIT_SUCCESS) {
      return return_with_mesg("Failed to run 'eval_psi_and_dpsidx_arr()");
    }
    

    //// Store `r_p_arr` and `v_p_arr`
    std::copy(
        r_p_arr_1d, r_p_arr_1d + r_p_arr_size,
        r_p_t_arr_1d + (i_t) * r_p_arr_size);
    std::copy(
        v_p_arr_1d, v_p_arr_1d + r_p_arr_size,
        v_p_t_arr_1d + (i_t) * r_p_arr_size);


    //// Update time index
    i_t++;

#else // BOHM
    
    // propagate one step forward in (real) time. 
    wf.propagate(timestep, time, g_prop, hamilton, me, staticpot, my_m_quantum_num, nuclear_charge);

#endif // BOHM
    

    if ((ts+1)%(ldumpwidth*10)==0) {
      cout << "timestep " << ts+1 << " of " << lno_of_ts+1 << ", Norm of wave function: " << N << endl;
    };
    if ((ts+1)% (wf_backup_period) == 0) {
      std::ofstream current_wf_bin_file(current_wf_bin_file_name, std::ios::binary);
      if (wf.dump_to_file_binary(current_wf_bin_file)) {
        std::cerr << "[ERROR] Failed to save wavefunction to " << current_wf_bin_file_name << endl;
        return 1; }
      current_wf_bin_file.close();
      cout << "[ LOG ] Wavefunction have been saved to " << current_wf_bin_file_name 
        << " at timestep " << ts+1 << " of " << lno_of_ts+1 << endl;
    }
  }; // end of real-time-propagation loop
#ifdef HAVE_BOOST
    cout << "time step took " << tim.elapsed() << " seconds" << endl;
#endif



#ifdef BOHM
  const size_t traj_file_size = N_t * r_p_arr_size * sizeof(double);
  string file_name_traj_r_t = string("traj-r-t.bin");
  string file_name_traj_v_t = string("traj-v-t.bin");

  std::ofstream file_traj_r_t;
  fprintf(stdout,"[ LOG ] Writing trajectory position to file: '%s' ... ", 
      file_name_traj_r_t.c_str());
  file_traj_r_t.open(file_name_traj_r_t, std::ios::binary);
  file_traj_r_t.write((char *) r_p_t_arr_1d, traj_file_size);
  file_traj_r_t.close();
  fprintf(stdout,"done\n");

  std::ofstream file_traj_v_t;
  fprintf(stdout,"[ LOG ] Writing trajectory velocity to file: '%s' ... ", 
      file_name_traj_v_t.c_str());
  file_traj_v_t.open(file_name_traj_v_t, std::ios::binary);
  file_traj_v_t.write((char *) v_p_t_arr_1d, traj_file_size);
  file_traj_v_t.close();
  fprintf(stdout,"done\n");

//  print_bar();
//
//  printf("r_p_t_arr: \n");
//  for (int i_t = 0; i_t < N_t; i_t++) {
//    print_3vec_p_arr(r_p_t_arr[i_t],N_p,N_r_dim);
//    printf("\n");
//  } 
//  print_bar();


  delete [] l_arr;
  delete [] m_arr;

  delete [] k_st_arr_1d;
  delete [] k_st_arr_2d;
  delete [] k_st_arr;
  delete [] r_p_i_st_arr_1d;
  delete [] r_p_i_st_arr;
  
  // Deallocate memory for arrays
  delete [] r_p_arr[0];
  delete [] r_p_arr;
  delete [] v_p_arr[0];
  delete [] v_p_arr;
//  delete [] psi_arr[0];
  delete [] psi_arr;
  delete [] r_p_t_arr_1d;
  delete [] r_p_t_arr;

#endif // BOHM



  fclose(file_obser);

  double yield_N = (1.0 - N);
  double yield_P = (1.0 - real(conj(P)*P));
  fprintf(file_yield, "%15.10le %15.10le\n", yield_N, yield_P);
  fclose(file_yield);

  wf.dump_to_file_sh(g_prop, file_wf, 1); // final wf is saved
  fclose(file_wf); 


  //// For split-step calculation
  std::ofstream current_wf_bin_file(current_wf_bin_file_name, std::ios::binary);
  wf.dump_to_file_binary(current_wf_bin_file);
  current_wf_bin_file.close();
  ////

 
  if (iv!=0) {
    fprintf(stdout, "%s is written.\n", str_fname_obser.c_str());
    fprintf(stdout, "%s is written.\n", str_fname_wf.c_str());
    if (start_time_index == 0) {
      fprintf(stdout, "%s is written.\n", str_fname_vpot.c_str()); }
    fprintf(stdout, "%s is written.\n", str_fname_logfi.c_str());
    fprintf(stdout, "%s is written.\n", str_fname_yield.c_str());
    fprintf(stdout, "Hasta la vista...\n");
    if (use_ppp) {
      cout << "[ LOG ] Run `PPP` for `time_surff`\n";
    }
  };
};
//
// end of main program
//
