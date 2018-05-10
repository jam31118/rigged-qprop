#include <iostream>
#include <complex>
#include <functional>
#include <string>

#include <bar.h>
#include <fluid.h>
#include <grid.h>
#include <hamop.h>
#include <wavefunction.h>
#include <powers.hh>
#include <parameter.hh>
#include <smallHelpers.hh>

// #include <../base/kbhit.cc> // kbhit() for unix 
//
// Functions determining the potentials are defined in potentials.hh
//

#include "potentials.hh"

using std::cout;
using std::endl;
using std::string;

void ell_m_consistency(long ell, long m, grid g) {
  if (ell<labs(m)) {
    cerr << "|m| is greater than ell" << endl;
    exit(-1);
  };
  if (g.ngps_y()<ell+1) {
    cerr << "grid too small for ell" << endl;
    exit(-1);
  };
};

bool nlm_are_consistent(long n, long ell, long m, grid g) {
  bool consistent = true;
  consistent &= n > ell;
  consistent &= ell >= labs(m);
  consistent &= g.ngps_y() > ell;
  consistent &= n >= 1; // for valid n
  return consistent;
}
  

int substract_component(grid g, wavefunction &wf1, wavefunction &wf2) {
  // DESCRIPTION //
  // Compute the wf2 component in wf1 
  // .. and substract it so that wf1 and wf2 are mutually orthogonal.
  // No normalization is done on the result.
  
  if (wf1.wf_size() != wf1.wf_size()) {
    fprintf(stderr, "[ERROR] The sizes of two input wavefunction are different.");
    return 1;
  }
  
  cplxd component_of_wf2_in_wf1 = wf1 * wf2 * g.delt_x();
  wf1 = wf1 - component_of_wf2_in_wf1 * wf2;

  return 0; 
}

// all wf in wf_arr should be initialized with size of g
int get_n_th_eigenstate(long n, const long my_ell_quantum_num, wavefunction *wf_arr, const cplxd timestep, 
    const long lno_of_ts, grid& g, hamop& hamilton, const int me, 
    wavefunction& staticpot, scalarpot& scalarpotx, const long my_m_quantum_num,
    FILE *file_obser_imag, int iv, double *acc_in, double *E_tot_in, const double acc_tol, int iinitmode) {

  cout << "[ LOG ] entered get-n-th with n: "<< n << endl;

  if ( n == my_ell_quantum_num+1 ) {} // just pass
  else if ( n > my_ell_quantum_num+1 ) {
    if (get_n_th_eigenstate(n-1, my_ell_quantum_num, wf_arr, timestep, lno_of_ts, g, hamilton, 
          me, staticpot, scalarpotx, my_m_quantum_num, file_obser_imag, iv, acc_in, E_tot_in, acc_tol, iinitmode) != 0) {
      std::cerr << "[ERROR] Failed to get " << n-1 << "-th eigenstate\n";
      return 1;
    }
  } else { 
    std::cerr << "[ERROR] Unexpected value for n: " << n << endl;
    return 1;
  }

  wavefunction *p_wf = &wf_arr[n-1];

//  int iinitmode = 2;
  fluid ell_init, m_init;
  ell_init[0] = my_ell_quantum_num;
  m_init[0] = my_m_quantum_num;
  (*p_wf).init(g.size());
  if (g.dimens()==34) { (*p_wf).init(g, iinitmode, 1.0, ell_init); } 
  else if (g.dimens() == 44) { (*p_wf).init_rlm(g, iinitmode, 1.0, ell_init, m_init); }
  else { 
    std::cerr << "[ERROR] Unexpected propagation mode: " << g.dimens() << endl;
    return 1;
  }
  (*p_wf).normalize(g);

  cout << "wf[0] = " << (*p_wf)[0] << endl;

  double E_tot_prev, acc, E_tot = *E_tot_in;
  fprintf(file_obser_imag,"[ LOG ] ****************************\n");
  fprintf(file_obser_imag,"[ LOG ] Constructing %ld-th eigenfucntion\n",n);

  if ((iv==1) && (is_time(0.1)))
    fprintf(stdout,
" time                   \
E_tot                  \
accuracy              \
step\n");

  for (long ts=0; ts<lno_of_ts; ts++) {
    const double time = double(ts)*imag(timestep);
//    cout << "[ LOG ] in loop, timestep: " << timestep << endl;
    // calculate the total energy
    E_tot_prev = E_tot;
    E_tot = real((*p_wf).energy(0.0, g, hamilton, me, staticpot, scalarpotx.get_nuclear_charge()));
    acc = fabs((E_tot_prev-E_tot)/(E_tot_prev+E_tot));
    fprintf(file_obser_imag,"% 20.15le % 20.15le %20.15le %ld\n", time, E_tot, acc, ts);

    if ( iv == 1 ) {
      fprintf(stdout,"% 20.15le % 20.15le %20.15le %ld\n", time, E_tot, acc, ts); }

//    cout << "[ LOG ] before propataiong at time: " << time << endl;
//    cout << "[ LOG ] g.dimens() = " << g.dimens() << endl;
//    cout << "[ LOG ] staticpot[1] = " << staticpot[1] << endl;
//    cout << "[ LOG ] my_m_quantum_num = " << my_m_quantum_num << endl;
//    cout << "[ LOG ] scalarpotx.get_nuclear_charge() = " << scalarpotx.get_nuclear_charge() << endl; 

    (*p_wf).propagate(timestep, 0.0, g, hamilton, me, staticpot, my_m_quantum_num, scalarpotx.get_nuclear_charge());
//    cout << "[ LOG ] after propataiong at time: " << time << endl;

    for (long wf_index=my_ell_quantum_num; wf_index<n-1; wf_index++) {
      substract_component(g, *p_wf, wf_arr[wf_index]);
    }
//    cout << "[ LOG ] after substraction at time: " << time << endl;
    (*p_wf).normalize(g);

    if ( acc < acc_tol ) { cout << "[ LOG ] acc tolerance reached. Terminating imag propagation. . .\n"; break; }
  }
  *acc_in = acc;
  *E_tot_in = E_tot;

  fprintf(file_obser_imag,"[ LOG ] ****************************\n");
  return 0;
}



int main(int argc, char **argv) {
  //
  // variables
  //
  double acc, E_tot_prev;

  grid g;
  hamop hamilton;
  wavefunction staticpot, E_i, wf;


  parameterListe para_ini("initial.param");
  
  
  double alpha;
  try { alpha = para_ini.getDouble("effpot-alpha"); }
  catch (std::exception&) { alpha = 0.0; } // default value - if alpha = 0.0 means no effective potential but just become coulumb potential.
  scalarpot scalarpotx(para_ini.getDouble("nuclear-charge"), para_ini.getDouble("pot-cutoff"), alpha);
// scalarpot scalarpotx(para_ini.getDouble("nuclear-charge"), para_ini.getDouble("pot-cutoff"));
  
  //
  // input
  //
  // *** declare the grid ***
  g.set_dim(para_ini.getLong("qprop-dim")); // 44 elliptical polariz., 34 linear polariz.
  const double delta_r = para_ini.getDouble("delta-r");
  double grid_radial_size = para_ini.getDouble("radial-grid-size");
  g.set_ngps(long(grid_radial_size/delta_r), para_ini.getLong("ell-grid-size"), 1);  // <--------------------------------- max. Anzahl in r-Richtung, in ell-Richtung, immer 1
  g.set_delt(delta_r);  // <-------------------------------- delta r
  g.set_offs(0, 0, 0);

  int iinitmode = 2; // 1 -- random, 2 -- hydrogenic wf.

  int iv        = 1; // verbosity of stdout
     
  //
  // prepare for propagation ...
  // 

  // Number of imaginary time steps
  const double acc_tol = 1.0e-14;  // [NOTE] Global config

  long max_num_of_timesteps_temp;
  try { max_num_of_timesteps_temp = para_ini.getLong("max-num-of-imag-timesteps"); }
  catch (std::exception&) { max_num_of_timesteps_temp = 200000; }
  const long lno_of_ts = max_num_of_timesteps_temp;  // [NOTE] Global config -- maximum number of timesteps
  fluid ell_init, m_init;
  ell_init.init(g.ngps_z());
  m_init.init(g.ngps_z());
  
  const long my_m_quantum_num=para_ini.getLong("initial-m");
  const long my_ell_quantum_num=para_ini.getLong("initial-ell");
  ell_m_consistency(my_ell_quantum_num, my_m_quantum_num, g);
  long initial_n;
  try { initial_n = para_ini.getLong("initial-n"); }
  catch (std::exception&) { initial_n = my_ell_quantum_num + 1; } // default value
  if ( initial_n < 1 ) { std::cerr << "[ERROR] 'initial-n' must be equal or bigger than 1\n"; }
  
  if ( ! nlm_are_consistent(initial_n, my_ell_quantum_num, my_m_quantum_num, g) ) {
    std::cerr << "[ERROR] nlm are not consistent\n";
    std::cerr << "[ LOG ] n: " << initial_n << " ell: " << my_ell_quantum_num << " m: " << my_m_quantum_num << endl;
    return 1;
  }

  ell_init[0] = my_ell_quantum_num; // populated l quantum number (needed for initialization only)  <---------------------------- 1s, 2p, 3d, ?????
  m_init[0] = my_m_quantum_num;     // populated m quantum number (needed for initialization only)
  const int me = 0; // dummy here

  E_i.init(g.ngps_z());

  // set the purly imaginary timestep
  const double imag_timestep=g.delt_x()/4.0;

  // the Hamiltonian
  imagpot imaginarypot(0, 0.0);
  hamilton.init(g, always_zero2, always_zero2, always_zero2, scalarpotx, always_zero5, always_zero5, imaginarypot, always_zero2);
  
  // this is the linear and constant part of the Hamiltonian
  staticpot.init(g.size()); 
  staticpot.calculate_staticpot(g, hamilton);

  // create some files with appropriate appendices
  string common_prefix("imag-prop");  // [NOTE] should be in global config
  string str_fname_logfi=common_prefix+string(".log");
  FILE* file_logfi = fopen_with_check(str_fname_logfi, "w");

  string str_fname_obser=common_prefix+string("-observ.dat");
  FILE* file_obser_imag = fopen_with_check(str_fname_obser,"w");

  string str_fname_wf_ini=common_prefix+string("-wf_ini.dat");
  FILE* file_wf_ini = fopen_with_check(str_fname_wf_ini,"w");

  string str_fname_wf_fin=common_prefix+string("-wf_fin.dat");
  FILE* file_wf_fin = fopen_with_check(str_fname_wf_fin,"w");

  if (iv!=0) {
    cout  <<  str_fname_logfi  << " will be (re)written." << endl
	  <<  str_fname_wf_ini  << " will be (re)written." << endl
	  <<  str_fname_obser  << " will be (re)written." << endl
	  <<  str_fname_wf_fin  << " will be (re)written." << endl;
  };
  
  // *** new initialization ***
  // *** the wavefunction array 
  wf.init(g.size()); 
  if (g.dimens()==34) { wf.init(g, iinitmode, 1.0, ell_init); } 
  else if (g.dimens() == 44) { wf.init_rlm(g, iinitmode, 1.0, ell_init, m_init); }
  else { std::cerr << "[ERROR] Unexpected propagation mode: " << g.dimens() << endl; }
  
  wf.normalize(g);

  fprintf(stdout, "Norm of KS-orbital: %le\n", wf.norm(g));
  
  //
  // write in log file
  //
  fprintf(file_logfi,"Imaginary-time propagation\n");
  fprintf(file_logfi,"Grid: \n");
  fprintf(file_logfi,"g.ngps_x() = %ld\n", g.ngps_x());
  fprintf(file_logfi,"g.ngps_y() = %ld\n", g.ngps_y());
  fprintf(file_logfi,"g.ngps_z() = %ld\n", g.ngps_z());
  fprintf(file_logfi,"g.dimens() = %d\n\n", g.dimens());
  fprintf(file_logfi,"g.delt_x() = %20.15le\n", g.delt_x());

  fprintf(file_logfi,"imag_timestep     = %20.15le\n", imag_timestep);
  fprintf(file_logfi,"lno_of_ts         = %ld\n", lno_of_ts);
  fprintf(file_logfi,"nuclear_charge    = %20.15le\n", scalarpotx.get_nuclear_charge());
  fprintf(file_logfi,"iinitmode         = %d\n", iinitmode);
  fprintf(file_logfi,"str_fname_wf_ini = %s\n", str_fname_wf_ini.c_str());
  fprintf(file_logfi,"str_fname_obser  = %s\n", str_fname_obser.c_str());
  fprintf(file_logfi,"str_fname_wf_fin = %s\n", str_fname_wf_fin.c_str());
  fflush(file_logfi);


  // *** the orbitals before relaxation
  wf.dump_to_file_sh(g, file_wf_ini, 1);
  fclose(file_wf_ini);

  // ********************************************************
  // ***** imaginary propagation ****************************
  // ********************************************************


  double E_tot = 0.0;
  const cplxd timestep(0.0, -1.0*imag_timestep);
  
  wavefunction *wf_arr = new wavefunction[initial_n];
  wavefunction *p_wf;
  for (long wf_index = 0; wf_index < initial_n; wf_index++) {
    p_wf = &wf_arr[wf_index];
    (*p_wf).init(g.size());
    if (g.dimens()==34) { (*p_wf).init(g, iinitmode, 1.0, ell_init); } 
    else if (g.dimens() == 44) { (*p_wf).init_rlm(g, iinitmode, 1.0, ell_init, m_init); }
    else { std::cerr << "[ERROR] Unexpected propagation mode: " << g.dimens() << endl; }
    (*p_wf).normalize(g);
  }

  cout << "[ LOG ] right before get-n-th ... \n";

  if (get_n_th_eigenstate(initial_n, my_ell_quantum_num, wf_arr, timestep, lno_of_ts, g, hamilton, 
          me, staticpot, scalarpotx, my_m_quantum_num, file_obser_imag, iv, &acc, &E_tot, acc_tol, iinitmode) != 0) {
    std::cerr << "[ERROR] Failed to get " << initial_n-1 << "-th eigenstate\n";
    return 1;
  }
  wf = wf_arr[initial_n-1];


//  // printf("Press any key to finish imaginary-time propagation...\n");
//  if ((iv==1) && (is_time(0.1)))
//    fprintf(stdout,
//" time                   \
//E_tot                  \
//accuracy              \
//step\n");
//
//  for (long ts=0; ts<lno_of_ts; ts++) {
//    const double time = double(ts)*imag(timestep);  
//    //
//    // calculate the total energy
//    // 
//    E_tot_prev = E_tot;
//    E_tot = real(wf.energy(0.0, g, hamilton, me, staticpot, scalarpotx.get_nuclear_charge()));
//    acc = fabs((E_tot_prev-E_tot)/(E_tot_prev+E_tot));
//    fprintf(file_obser_imag,"% 20.15le % 20.15le %20.15le %ld\n", time, E_tot, acc, ts);
//
//    if ((iv==1)) {
//      fprintf(stdout,"% 20.15le % 20.15le %20.15le %ld\n", time, E_tot, acc, ts); }
//
//    wf.propagate(timestep, 0.0, g, hamilton, me, staticpot, my_m_quantum_num, scalarpotx.get_nuclear_charge());
//
//    wf.normalize(g);
//
//  };
 
  fprintf(file_logfi, "acc = %le\n", acc);

  fclose(file_obser_imag);
  fclose(file_logfi);

  wf.dump_to_file_sh(g,file_wf_fin,1);
  fclose(file_wf_fin);


  //// Added for initial energy output
  string filename_initial_energy = string("initial-energy.dat");
  FILE *file_initial_energy = fopen_with_check(filename_initial_energy, "w");
  fprintf(file_initial_energy, "%20.15le", E_tot);
  fclose(file_initial_energy);
  //// END


  //// Added for regriding to propagation grid ////
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");
  // *** declare the grid for propagation ***
  grid g_prop;
  g_prop.set_dim(para_ini.getLong("qprop-dim"));
  double beyond_R_distance_temp;
  try { beyond_R_distance_temp = para_tsurff.getDouble("beyond-R"); }
  catch (std::exception&) { beyond_R_distance_temp = 0.0; }
  const double beyond_R_distance = beyond_R_distance_temp;
  double grid_size = para_tsurff.getDouble("R-tsurff") + beyond_R_distance + para_prop.getDouble("imag-width");

  // [NOTE] It assumes there's no quiver-amplitude addition.
  // However, it should become possible to add more space beyond R-tsurff
  // Let's make a parameter for this.
//  bool add_quiver_ampl;
//  const double omega = para_prop.getDouble("omega");
//  double E_0 = para_prop.getDouble("max-electric-field");
//  const double quiver_amplitude = E_0/pow2(omega);
//  try { add_quiver_ampl = para_prop.getBool("add-quiver-amplitude-in-real-prop"); }
//  catch (std::exception& e) { add_quiver_ampl = false; }
//  cout << "add_quiver_ampl: " << add_quiver_ampl << endl;
//  if (add_quiver_ampl) { grid_size += quiver_amplitude; }
  //
  g_prop.set_ngps(long(grid_size/delta_r), para_prop.getLong("ell-grid-size"), 1); 
  g_prop.set_delt(delta_r);
  g_prop.set_offs(0, 0, 0);

  wavefunction wf_prop; 
  wf_prop.init(g_prop.size()); 
  wf_prop.regrid(g_prop, g, wf);

  //// END ////


  //// 180426 Added for split-time calc
  string current_wf_bin_file_name = string("ini-wf.bin");
  std::ofstream current_wf_bin_file(current_wf_bin_file_name, std::ios::binary);
  wf_prop.dump_to_file_binary(current_wf_bin_file);
  current_wf_bin_file.close();
  ////


  if (iv!=0) {
  cout  <<  str_fname_logfi  << " is written." << endl
	<<  str_fname_wf_ini  << " is written." << endl
	<<  str_fname_obser  << " is written." << endl
	<<  str_fname_wf_fin  << " is written." << endl;
  };

  fprintf(stdout, "Hasta la vista...\n");
};

// 
// end of main program
// 
