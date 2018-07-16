#include "real-prop.h"

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


  // How long do the slowest electrons have time to reach the t-SURFF boundary
  // [NOTE] Keep in mind that the electron doesn't have a position 
  // and there's no concept like 'reach somewhere' for the electron especially in this quantum treatment.
  const double time_surff=para_tsurff.getDouble("R-tsurff")/para_tsurff.getDouble("p-min-tsurff");
  const double duration=pulse_duration+time_surff;
  cout << "[ LOG ] pulse_duration : " << pulse_duration << endl;
  cout << "[ LOG ] time_surff : " << time_surff << endl;
  
  
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
  string initial_wf_file_name = string("ini-wf.bin"); // [NOTE] global variable
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
  if (start_time_index == 0) { start_wf_file_name = initial_wf_file_name; } // [NOTE] global variable
  else if (start_time_index > 0) { start_wf_file_name = string("current-wf.bin"); } // [NOTE] global variable
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
  for (long ts=start_time_index; ts < max_time_index; ts++) {
    const double time=real_timestep*double(ts);
    // save the orbitals \varphi_{\ell}(\RI) and the derivative \partial_r\varphi_{\ell}(r)|_{r=\RI}
    tsurff_save_wf(wf);

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
    //
    // propagate one step forward in (real) time.
    // 
    wf.propagate(timestep, time, g_prop, hamilton, me, staticpot, my_m_quantum_num, nuclear_charge);

    if (ts%(ldumpwidth*10)==0) {
      cout << "timestep " << ts << " of " << lno_of_ts << ", Norm of wave function: " << N << endl;
    };
    if (ts % (wf_backup_period) == 0) {
      std::ofstream current_wf_bin_file(current_wf_bin_file_name, std::ios::binary);
      if (wf.dump_to_file_binary(current_wf_bin_file)) {
        std::cerr << "[ERROR] Failed to save wavefunction to " << current_wf_bin_file_name << endl;
        return 1; }
      current_wf_bin_file.close();
      cout << "[ LOG ] Wavefunction have been saved to " << current_wf_bin_file_name 
        << " at timestep " << ts << " of " << lno_of_ts << endl;
    }
  }; // end of real-time-propagation loop
#ifdef HAVE_BOOST
    cout << "time step took " << tim.elapsed() << " seconds" << endl;
#endif

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
  };
};
//
// end of main program
//
