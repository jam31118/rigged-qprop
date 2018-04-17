#include <iostream>
#include <memory>
#include <complex>
typedef std::complex<double> cplxd;
typedef std::unique_ptr<cplxd[]> cplxd_ptr;

#include <grid.h>
#include <hamop.h>
#include <wavefunction.h>
#include <parameter.hh>
#ifdef HAVE_BOOST
#include <boost/timer.hpp>
#endif
#include <smallHelpers.hh>
#include <powers.hh>
#include <tsurffSpectrum.hh>

// Functions, which determine potentials
#include "potentials.hh"

using std::endl;
using std::cout;

void print_banner() {
  fprintf(stdout, " IONIZATION: hydrogen\n");
  fprintf(stdout, " (C) Copyright by Bauer D and Koval P, Heidelberg (2005)\n");
  fprintf(stdout, " -------------------------------------------------------\n");
};

int main(int argc, char **argv) {
  grid g_prop, g_load;
  wavefunction staticpot, wf, wf_load;

  print_banner();
  
  // input (result from imaginary propagation)
  string str_fname_wf_ini=string("./hydrogen_im-0-wf_fin.dat");
  FILE* file_wf_ini = fopen_with_check(str_fname_wf_ini, "r");
  
  int isave_wf      = 1; // 1 -- save or 0 -- don't save wavefunctions
  int iv            = 1; // verbosity of stdout

  // get parameters
  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");
  
  // *** declare the grid for load ***
  g_load.set_dim(para_ini.getLong("qprop-dim"));
  const double delta_r = para_ini.getDouble("delta-r");
  g_load.set_ngps(long(para_ini.getDouble("radial-grid-size")/delta_r), para_ini.getLong("ell-grid-size"), 1);
  g_load.set_delt(delta_r, 0.0, 0.0);
  g_load.set_offs(0, 0, 0);
  
  const int my_m_quantum_num = para_ini.getLong("initial-m");
  const double nuclear_charge = para_ini.getDouble("nuclear-charge");

  // everything related to the laser
  const double omega = para_prop.getDouble("omega");
  double E_0     = para_prop.getDouble("max-electric-field");
  const double quiver_amplitude = E_0/pow2(omega);
  double I_p     = nuclear_charge*nuclear_charge / 2.0; // warn: hydrogenic energy!
  // only A_z in this example
  vecpot vecpot_x(omega, 1.0, 0.0, 0.0);
  vecpot vecpot_y(omega, 1.0, 0.0, 0.0);
  vecpot vecpot_z(omega, para_prop.getDouble("num-cycles"), E_0, 0.0);
  const double U_p     = vecpot_z.get_Up();
  const double gamma_K = sqrt(I_p / 2.0 / U_p);     
  const double pulse_duration=vecpot_z.get_duration();
  // how long do the slowest electrons have time to reach the t-SURFF boundary
  const double time_surff=para_tsurff.getDouble("R-tsurff")/para_tsurff.getDouble("p-min-tsurff");
  const double duration=pulse_duration+time_surff;
  
  // *** declare the grid for propagation ***
  g_prop.set_dim(para_ini.getLong("qprop-dim"));
  const double grid_size=para_prop.getDouble("imag-width")+para_tsurff.getDouble("R-tsurff")+quiver_amplitude;
  g_prop.set_ngps(long(grid_size/delta_r), para_prop.getLong("ell-grid-size"), 1); 
  g_prop.set_delt(delta_r);
  g_prop.set_offs(0, 0, 0);
  
  // output that will be created by hydrogen_re
  string common_prefix("hydrogen_re");
  string str_fname_logfi=common_prefix+string(".log");
  FILE* file_logfi = fopen_with_check(str_fname_logfi, "w");
  string str_fname_yield=common_prefix+string("-yield.dat");
  FILE* file_yield = fopen_with_check(str_fname_yield, "w");
  string str_fname_obser=common_prefix+string("-obser.dat");
  FILE* file_obser = fopen_with_check(str_fname_obser, "w");
  
  if (iv!=0) {
    fprintf(stdout, "%s will be (re)written.\n", str_fname_yield.c_str());
    fprintf(stdout, "%s will be (re)written.\n", str_fname_logfi.c_str());
  };

  // create an instance of the class for doing the tsurff related work
  tsurffSaveWF tsurff_save_wf(para_ini, para_prop, para_tsurff, g_prop);

  // the absorbing imaginary potential
  const long imag_potential_width=long(para_prop.getDouble("imag-width")/delta_r);
  imagpot imaginarypot(imag_potential_width);
  // set the binding potential and the hamiltonian
  scalarpot scalarpotx(nuclear_charge, para_ini.getDouble("pot-cutoff"));
  hamop hamilton;
  hamilton.init(g_prop, vecpot_x, vecpot_y, vecpot_z, scalarpotx, always_zero5, always_zero5, imaginarypot, always_zero2);

  // this is the linear and constant part of the Hamiltonian
  staticpot.init(g_prop.size()); 
  staticpot.calculate_staticpot(g_prop, hamilton);

  // *** wavefunction array 
  wf.init(g_prop.size()); 
  wf_load.init(g_load.size());
  
  // *** wavefunction initialization ***
  wf_load.init(g_load, file_wf_ini, 0, iv);
  wf.regrid(g_prop, g_load, wf_load);    
  fclose(file_wf_ini);

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
  fprintf(file_logfi, "str_fname_wf_ini = %s\n", str_fname_wf_ini.c_str());
  fprintf(file_logfi, "str_fname_obser  = %s\n", str_fname_obser.c_str());
  if (isave_wf==1) 
    fprintf(file_logfi, "str_fname_wf = %s\n", str_fname_wf.c_str());
  // fprintf(file_logfi, "n_c        = %lf\n", n_c);  
  // fprintf(file_logfi, "n_r        = %lf\n", n_r);
  fprintf(file_logfi, "omega        = %15.10le\n", omega);
  // fprintf(file_logfi, "phi_cep      = %15.10le\n", phi_cep);
  fprintf(file_logfi, "duration = %15.10le\n", duration);
  // fprintf(file_logfi, "I_0      = %15.10le, W/cm^2\n", lintens);
  fprintf(file_logfi, "E_0      = %15.10le\n", E_0);
  fprintf(file_logfi, "gamma_K  = %15.10le\n", gamma_K);
  fprintf(file_logfi, "gamma_K      = %15.10le ", gamma_K);
  if (gamma_K>1.0) {
    fprintf(file_logfi, "(multi-photon regime)\n");
  }
  else {
    fprintf(file_logfi, "(tunneling regime)\n");
  };
  fflush(file_logfi);
  fclose(file_logfi);

  long ldumpwidth(1.0/real_timestep);  // output once every a.u. of time
  int me = 0; // dummy here
  // write vector potential to file
  string str_fname_vpot_z=common_prefix+string("-vpot_z.dat");
  FILE* file_vpot_z = fopen_with_check(str_fname_vpot_z, "w");
  for (long ts=0; ts<lno_of_ts; ts++) {
    const double time=real_timestep*double(ts);
    if (ts%ldumpwidth==0)
      fprintf(file_vpot_z, "%15.10le %15.10le\n", time, vecpot_z(time, me));
  };
  fclose(file_vpot_z);

  // ********************************************************
  // ***** real time propagation ****************************
  // ********************************************************
  cplxd timestep=cplxd(real_timestep, 0.0);
  cplxd P;
  double N;
#ifdef HAVE_BOOST
  boost::timer tim;
#endif
  for (long ts=0; ts<lno_of_ts; ts++) {
#ifdef HAVE_BOOST
    tim.restart();
#endif
    const double time=real_timestep*double(ts);
    // save the orbitals \varphi_{\ell}(\RI) and the derivative \partial_r\varphi_{\ell}(r)|_{r=\RI}
    tsurff_save_wf(wf);

    if (ts%ldumpwidth==0) {
      // calculate total energy, projection onto initial state, norm, and <z>
      double E_tot = real(wf.energy(0.0, g_prop, hamilton, me, staticpot, nuclear_charge));
      P = wf.project(g_prop, g_load, wf_load, 0);
      N = wf.norm(g_prop);
      double z_expect = real(wf.expect_z(g_prop));
      fprintf(file_obser, "%15.10le %15.10le %15.10le %15.10le %15.10le\n", time, E_tot, real(conj(P)*P), N, z_expect);
    };
    //
    // propagate one step forward in (real) time.
    // 
    wf.propagate(timestep, time, g_prop, hamilton, me, staticpot, my_m_quantum_num, nuclear_charge);

    if (ts%(ldumpwidth*10)==0) {
      cout << "timestep " << ts << " of " << lno_of_ts << ", Norm of wave function: " << N << endl;
    };
#ifdef HAVE_BOOST
    cout << "time step took " << tim.elapsed() << " seconds" << endl;
#endif
  }; // end of real-time-propagation loop

  fclose(file_obser);

  double yield_N = (1.0 - N);
  double yield_P = (1.0 - real(conj(P)*P));
  fprintf(file_yield, "%15.10le %15.10le\n", yield_N, yield_P);
  fclose(file_yield);

  wf.dump_to_file_sh(g_prop, file_wf, 1); // final wf is saved
  fclose(file_wf); 
 
  if (iv!=0) {
    fprintf(stdout, "%s was read.\n",   str_fname_wf_ini.c_str());
    fprintf(stdout, "%s is written.\n", str_fname_obser.c_str());
    fprintf(stdout, "%s is written.\n", str_fname_wf.c_str());
    fprintf(stdout, "%s is written.\n", str_fname_vpot_z.c_str());
    fprintf(stdout, "%s is written.\n", str_fname_logfi.c_str());
    fprintf(stdout, "%s is written.\n", str_fname_yield.c_str());
    fprintf(stdout, "Hasta la vista...\n");
  };
};
//
// end of main program
//
