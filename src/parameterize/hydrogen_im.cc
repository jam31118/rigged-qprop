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
  

int main(int argc, char **argv) {
  //
  // variables
  //
  double acc, E_tot_prev;

  grid g;
  hamop hamilton;
  wavefunction staticpot, E_i, wf;

  fprintf(stdout," HYDROGEN_IM: Imaginary-time propagation for hydrogen atom\n");
  fprintf(stdout," (C) Copyright by Bauer D and Koval P, Heidelberg (2005)\n");
  fprintf(stdout," -------------------------------------------------------\n");

  parameterListe para_ini("initial.param");
  
  scalarpot scalarpotx(para_ini.getDouble("nuclear-charge"), para_ini.getDouble("pot-cutoff"));
  
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
  const long lno_of_ts = 5000;  // [NOTE] Global config
  fluid ell_init, m_init;
  ell_init.init(g.ngps_z());
  m_init.init(g.ngps_z());
  
  const long my_m_quantum_num=para_ini.getLong("initial-m");
  const long my_ell_quantum_num=para_ini.getLong("initial-ell");
  ell_m_consistency(my_ell_quantum_num, my_m_quantum_num, g);
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
  string common_prefix("hydrogen_im-");
  string str_fname_logfi=common_prefix+to_string(my_m_quantum_num)+string(".log");
  FILE* file_logfi = fopen_with_check(str_fname_logfi, "w");

  string str_fname_obser=common_prefix+to_string(my_m_quantum_num)+string("-observ.dat");
  FILE* file_obser_imag = fopen_with_check(str_fname_obser,"w");

  string str_fname_wf_ini=common_prefix+to_string(my_m_quantum_num)+string("-wf_ini.dat");
  FILE* file_wf_ini = fopen_with_check(str_fname_wf_ini,"w");

  string str_fname_wf_fin=common_prefix+to_string(my_m_quantum_num)+string("-wf_fin.dat");
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

  // printf("Press any key to finish imaginary-time propagation...\n");
  if ((iv==1) && (is_time(0.1)))
    fprintf(stdout,
" time                   \
E_tot                  \
accuracy              \
step\n");

  for (long ts=0; ts<lno_of_ts; ts++) {
    const double time = double(ts)*imag(timestep);  
    //
    // calculate the total energy
    // 
    E_tot_prev = E_tot;
    E_tot = real(wf.energy(0.0, g, hamilton, me, staticpot, scalarpotx.get_nuclear_charge()));
    acc = fabs((E_tot_prev-E_tot)/(E_tot_prev+E_tot));
    fprintf(file_obser_imag,"% 20.15le % 20.15le %20.15le %ld\n", time, E_tot, acc, ts);

    if ((iv==1)) {
      fprintf(stdout,"% 20.15le % 20.15le %20.15le %ld\n", time, E_tot, acc, ts); }

    wf.propagate(timestep, 0.0, g, hamilton, me, staticpot, my_m_quantum_num, scalarpotx.get_nuclear_charge());

    wf.normalize(g);

  };
 
  fprintf(file_logfi, "acc = %le\n", acc);

  fclose(file_obser_imag);
  fclose(file_logfi);

  wf.dump_to_file_sh(g,file_wf_fin,1);
  fclose(file_wf_fin);


  //// Added for regriding to propagation grid ////
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");
  // *** declare the grid for propagation ***
  grid g_prop;
  g_prop.set_dim(para_ini.getLong("qprop-dim"));
  double grid_size = para_prop.getDouble("imag-width") + para_tsurff.getDouble("R-tsurff");

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
