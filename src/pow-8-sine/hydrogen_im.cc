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

  grid g;
  hamop hamilton;
  wavefunction staticpot, E_i, wf;

  fprintf(stdout," HYDROGEN_IM: Imaginary-time propagation for hydrogen atom\n");
  fprintf(stdout," (C) Copyright by Bauer D and Koval P, Heidelberg (2005)\n");
  fprintf(stdout," -------------------------------------------------------\n");

  parameterListe para("initial.param");
  
  scalarpot scalarpotx(para.getDouble("nuclear-charge"), para.getDouble("pot-cutoff"));
  
  //
  // input
  //
  // *** declare the grid ***
  g.set_dim(para.getLong("qprop-dim")); // 44 elliptical polariz., 34 linear polariz.
  const double delta_r = para.getDouble("delta-r");
  g.set_ngps(long(para.getDouble("radial-grid-size")/delta_r), para.getLong("ell-grid-size"), 1);  // <--------------------------------- max. Anzahl in r-Richtung, in ell-Richtung, immer 1
  g.set_delt(delta_r);  // <-------------------------------- delta r
  g.set_offs(0, 0, 0);

  int iinitmode = 2; // 1 -- random, 2 -- hydrogenic wf.

  int iv        = 1; // verbosity of stdout
     
  //
  // prepare for propagation ...
  // 

  // Number of imaginary time steps
  const long lno_of_ts = para.getLong("max-imag-ts");
  const double wf_convergence_crit = para.getDouble("wf-convergence-crit");
  
  fluid ell_init, m_init;
  ell_init.init(g.ngps_z());
  m_init.init(g.ngps_z());
  
  const long my_m_quantum_num=para.getLong("initial-m");
  const long my_ell_quantum_num=para.getLong("initial-ell");
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
  string common_prefix("hydrogen_im");
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
  if (g.dimens()==34) 
    wf.init(g, iinitmode, 1.0, ell_init);
  else
    wf.init_rlm(g, iinitmode, 1.0, ell_init, m_init); // appropriate for prop. mode 44
  
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
  double acc(1.0), E_tot_prev, convergence_wf(1.0);
  wavefunction old_wf=wf;

  //  if ((iv==1) && (is_time(0.1)))
  if ((iv==1))
    fprintf(stdout,
" time                   \
E_tot                 \
energy convergence    \
eigenstate convergence \
step\n");

  for (long ts=0; ts<lno_of_ts; ts++) {
    const double time = double(ts)*imag(timestep);
    // propagate the wavefunction one imaginary timestep
    wf.propagate(timestep, 0.0, g, hamilton, me, staticpot, my_m_quantum_num, scalarpotx.get_nuclear_charge());
    // renormalize it
    wf.normalize(g);
    // compare to result of last step
    old_wf=old_wf-wf;
    convergence_wf=old_wf.norm(g);
    old_wf=wf;
    
    // calculate the total energy
    E_tot = real(wf.energy(0.0, g, hamilton, me, staticpot, scalarpotx.get_nuclear_charge()));
    acc = fabs((E_tot_prev-E_tot)/(E_tot_prev+E_tot));
    E_tot_prev=E_tot;
    fprintf(file_obser_imag, "% 20.15le % 20.15le %20.15le %20.15le %ld\n", time, E_tot, acc, convergence_wf, ts);

    //    if ((iv==1) && (is_time(0.1)))
    if ((iv==1))
      fprintf(stdout,"% 20.15le % 20.15le %20.15le %20.15le %ld\n", time, E_tot, acc, convergence_wf, ts);
    // stop if wave function has converged
    if (fabs(convergence_wf)<fabs(wf_convergence_crit) && ts>1e3) break;
  };
 
  fprintf(file_logfi, "energy convergence = %le\n", acc);
  fprintf(file_logfi, "eigenstate convergence = %le\n", convergence_wf);
  
  fclose(file_obser_imag);
  fclose(file_logfi);

  wf.dump_to_file_sh(g, file_wf_fin, 1);
  fclose(file_wf_fin);

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
