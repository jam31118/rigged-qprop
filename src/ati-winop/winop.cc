#include <iostream>
#include <vector>
#include <complex>
typedef std::complex<double> cplxd;

#include <wavefunction.h>
#include <fluid.h>
#include <grid.h>
#include <hamop.h>
#include <winop.h>
#include <smallHelpers.hh>
#include <gsl/gsl_sf_legendre.h>
#include <parameter.hh>
#include <powers.hh>
// definitions of potentials
#include "potentials.hh"

using std::vector;
using std::endl;
using std::cout;

void print_banner() {
  fprintf(stdout, " WINOP: Calculation of one-electron spectra\n");
  fprintf(stdout, " (C) Copyright by Bauer D and Koval P, Heidelberg (2005)\n");
  fprintf(stdout, " -------------------------------------------------------\n");
};

int main(int argc, char **argv) {
  print_banner();
  // set verbosity
  const int iv=1;

  // get parameters
  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_winop("winop.param");
  
  //
  // set a few essential parameters for the calculation of the spectra
  //
  const int lnoofwfoutput=1;      // 0 -- initial state, 1 -- final state.
  const long num_energy=max(1l, para_winop.getLong("num-energy"));   // number of energy bins
  const double energy_max=para_winop.getDouble("energy-max"); // half width of the window operator (72)
  const double energy_min=para_winop.getDouble("energy-min");
  const double delta_energy=(energy_max-energy_min)/double(num_energy);
  const double gamma=0.5*delta_energy; // width of the bin is 2*gamma
  cout << "gamma: " << gamma << endl;
  const double w_norm=7.0*M_PI*gamma/(32.0*sin(M_PI/8.0));
  // define for how many angles the spectrum should be evaluated 
  const long num_theta=max(1l, para_winop.getLong("num-theta"));
  const long num_phi=max(1l, para_winop.getLong("num-phi"));
  const long num_dir=num_theta*num_phi;

  const int my_m_quantum_num=para_ini.getLong("initial-m");

  //
  // specify where the wavefunction (to be analyzed) is located
  //
  string str_fname_wf=string("hydrogen_re-wf.dat");
  FILE* file_wf=fopen_with_check(str_fname_wf, "r");

  // specific angles
  vector<double> theta(num_theta, 0.0);
  vector<double> cos_theta(num_theta, 0.0);
  for(long ltheta=0; ltheta<num_theta; ltheta++) {
    theta[ltheta]=(num_theta==1)?0.0:ltheta*M_PI/(num_theta-1.0);
    cos_theta[ltheta]=cos(theta[ltheta]);
  };
  
  vector<double> phi(num_phi, 0.0);
  for(long lphi=0; lphi<num_phi; lphi++) 
    phi[lphi]=lphi*2*M_PI/num_phi;

  const double delta_r=para_ini.getDouble("delta-r");
  const long qprop_dim=para_ini.getLong("qprop-dim");
  const long ell_grid_size=para_prop.getLong("ell-grid-size");

  grid  g_angle, g_load, g;
  // this grid is used to store the Y_lm; do not touch this !!!
  g_angle.set_dim(44);  // 44 works for both propagation mode 34 and 44
  g_angle.set_ngps(num_dir, ell_grid_size, 1);
  g_angle.set_delt(0.0);
  g_angle.set_offs(0, 0, 0);

  // this is the grid for the wavefunction that is to be analyzed
  const double R_max=para_prop.getDouble("R-max");
  const double grid_size=para_prop.getDouble("imag-width")+R_max;
  g_load.set_dim(qprop_dim);
  g_load.set_ngps(long(grid_size/delta_r), ell_grid_size, 1);
  g_load.set_delt(delta_r);
  g_load.set_offs(0, 0, 0);

  // this is the grid on which the wavefunction is analyzed
  // delt_r must be the same for g and g_load
  // increasing ngps_x improves the resolution in the continuum 
  g.set_dim(qprop_dim);
  g.set_ngps(para_winop.getLong("winop-radial-grid-size"), ell_grid_size, 1); 
  g.set_delt(delta_r);
  g.set_offs(0, 0, 0);

  wavefunction fullchi, result_lsub;
  if (qprop_dim==34) {
    fullchi.init(g.ngps_x()*g.ngps_y()); 
    result_lsub.init(g.ngps_y());
  }
  else if (qprop_dim==44) {
    fullchi.init(g.ngps_x()*g.ngps_y()*g.ngps_y()); 
    result_lsub.init(g.ngps_y()*g.ngps_y());
  }
  else { 
    fprintf(stderr, "err: wrong qprop-dim!=34 and qprop-dim!=44!\n"); exit(12);
  };

  fluid  V_ee_0;
  V_ee_0.init(g.ngps_x());

  //
  // the Hamiltonian
  //
  const double nuclear_charge=para_ini.getDouble("nuclear-charge");
  scalarpot scalarpotx(nuclear_charge, para_ini.getDouble("pot-cutoff"));
  hamop hamilton;
  hamilton.init(g, always_zero2, always_zero2, always_zero2, scalarpotx, always_zero5, always_zero5, always_zero_imag, always_zero2);

  // *** the wavefunction arrays
  wavefunction  wf, wf_load;
  wf.init(g.size());
  wf_load.init(g_load.size());

  // ***** calculate all the ylm's needed
  wavefunction ylm_array;
  ylm_array.init(g_angle.size());
  // Function: double gsl_sf_legendre_sphPlm (int l, int m, double x)
  // These routines compute the normalized associated Legendre polynomial \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x) suitable for use in spherical harmonics. The parameters must satisfy m >= 0, l >= m, |x| <= 1. Theses routines avoid the overflows that occur for the standard normalization of P_l^m(x)
  if (g.dimens()==34) {
    for (long l_index=my_m_quantum_num; l_index<g_angle.ngps_y(); l_index++) {      
      for (long ltheta=0; ltheta<num_theta; ltheta++) {
	for (long lphi=0; lphi<num_phi; lphi++) {
	  const long ylm_index=g_angle.rlmindex(ltheta+lphi*num_theta, l_index, 0);
	  if (my_m_quantum_num>=0)
	    ylm_array[ylm_index]=gsl_sf_legendre_sphPlm (l_index, my_m_quantum_num, cos_theta[ltheta])*exp(cplxd(0.0, double(my_m_quantum_num)*phi[lphi]));
	  else
	    ylm_array[ylm_index]=((my_m_quantum_num%2==0)?1.0:-1.0)*gsl_sf_legendre_sphPlm (l_index, labs(my_m_quantum_num), cos_theta[ltheta])*exp(cplxd(0.0, double(my_m_quantum_num)*phi[lphi]));
	};
      };
    };
  }
  else {
    for (long l_index=0; l_index<g_angle.ngps_y(); l_index++) {
      const long m_limit=l_index;
      for (long m_index=-m_limit; m_index<=m_limit; m_index++) {
	for (long ltheta=0; ltheta<num_theta; ltheta++) {
	  for (long lphi=0; lphi<num_phi; lphi++) {
	    const long ylm_index=g_angle.rlmindex(ltheta+lphi*num_theta, l_index, m_index);
	    if (m_index>=0)
	      ylm_array[ylm_index]=gsl_sf_legendre_sphPlm (l_index, m_index, cos_theta[ltheta])*exp(cplxd(0.0, double(m_index)*phi[lphi]));
	    else
	      ylm_array[ylm_index]=((m_index%2==0)?1.0:-1.0)*gsl_sf_legendre_sphPlm (l_index, labs(m_index), cos_theta[ltheta])*exp(cplxd(0.0, double(m_index)*phi[lphi]));
	  };
	};
      };
    };
  };
  
  if (iv==1) cout << "Reading ... " << endl;

  wf_load.init(g_load, file_wf, lnoofwfoutput, iv);
  fclose(file_wf);

  if (iv==1) cout << "... done." << endl;
  
  wf.regrid(g, g_load, wf_load);
  
  // Stdout
  if (iv==1) {
    fprintf(stdout, "Grid: \n");
    fprintf(stdout, "g.ngps_x()=%ld\n", g.ngps_x());
    fprintf(stdout, "g.ngps_y()=%ld\n", g.ngps_y());
    fprintf(stdout, "g.ngps_z()=%ld\n", g.ngps_z());
    fprintf(stdout, "g.dimens()=%d\n", g.dimens());
    fprintf(stdout, "g.delt_x()=%20.15le\n", g.delt_x());
    fprintf(stdout, "nuclear_charge   =%20.15le\n", nuclear_charge);
    fprintf(stdout, "str_fname_wf=%s\n", str_fname_wf.c_str());
    fflush(stdout);
  };

  wavefunction staticpot;
  staticpot.init(g.size());
  staticpot.calculate_staticpot(g, hamilton);
  
  // output files
  ofstream file_res(string("spectrum_")+to_string(my_m_quantum_num)+string(".dat"));
  file_res.precision(15);
  ofstream file_spectrum_polar(string("spectrum_polar")+to_string(my_m_quantum_num)+string(".dat"));
  file_spectrum_polar.precision(17);

  // *********************************
  // *** run over the energies
  // ********************************* 
  for (long lenergy=0; lenergy<num_energy; lenergy++) {
    const double energy=energy_min+lenergy*delta_energy;
    file_res << energy << " " << sqrt(2.0*energy) << " " ;

    cplxd result_tot(0.0, 0.0);
    winop_fullchi(fullchi, result_lsub, &result_tot, energy, gamma, staticpot, V_ee_0, nuclear_charge, g, wf, iv);

    // write the partial results (i.e., for individual l (and m)) to file
    for(long l_index=0; l_index<g.ngps_y(); l_index++) {
      if (g.dimens()==44) { 
	const long m_limit=l_index;
	for(long m_index=-m_limit; m_index<=m_limit; m_index++) {
	  file_res << real(result_lsub[((l_index+1)*(l_index+1)-(l_index+1)+m_index)])/w_norm << " ";
	};
      }
      else { // linear polarization and fixed value of m
	file_res << real(result_lsub[l_index])/w_norm << " ";
      };
    };
    // write the total result for the energy resolved spectrum to file
    file_res << real(result_tot)/w_norm << " ";

    // space for the angle-dependent result
    wavefunction angle_result;
    angle_result.init(num_dir);      
    angle_result.nullify();
    
    // calculate angle-dependent spectrum (equation (77) in original QPROP paper)
    if (g.dimens()==34) { // linear polarization and fixed value of m
      const long lphi=0;
      for (long ltheta=0; ltheta<num_theta; ltheta++) {
	cplxd summingup(0.0, 0.0);
	// loop for radial integration
	for (long rindex=0; rindex<g.ngps_x(); rindex++) {
	  summingup=cplxd(0.0, 0.0);
	  for (long l_index=my_m_quantum_num; l_index<g.ngps_y(); l_index++) {
	    const long index=g.index(rindex, l_index, 0, 0);
	    const cplxd b_El=fullchi[index]*ylm_array[g_angle.rlmindex(ltheta+num_theta*lphi, l_index, 0)];
	    summingup+=b_El;
	  };
	  angle_result[ltheta+num_theta*lphi]+=summingup*conj(summingup);
	};
	angle_result[ltheta+num_theta*lphi]*=g.delt_x();
	// write result of integration for every value of \theta
	file_spectrum_polar << energy << " " << sqrt(2.0*energy) << " " << theta[ltheta] << " " << real(angle_result[ltheta+num_theta*lphi])/w_norm << endl;
      };
    }
    else { // polarization in xy-plane
      for (long ltheta=0; ltheta<num_theta; ltheta++) {
	for (long lphi=0; lphi<num_phi; lphi++) {
	  cplxd summingup(0.0, 0.0);
	  // loop for radial integration
	  for (long rindex=0; rindex<g.ngps_x(); rindex++) {
	    summingup=cplxd(0.0,0.0);
	    for (long l_index=0; l_index<g.ngps_y(); l_index++) {
	      const long m_limit=l_index;
	      for (long m_index=-m_limit; m_index<=m_limit; m_index++) {
		const long index=g.index(rindex, l_index, m_index, 0);
		summingup+=fullchi[index]*ylm_array[g_angle.rlmindex(ltheta+num_theta*lphi, l_index, m_index)];
	      };
	    };
	    angle_result[ltheta+num_theta*lphi]+=summingup*conj(summingup);
	  };
	  angle_result[ltheta+num_theta*lphi]*=g.delt_x();
	  // write result of integration for every value of \varphi
	file_spectrum_polar << energy << " " << sqrt(2.0*energy) << " " << theta[ltheta] << " " << phi[lphi] << " " << real(angle_result[ltheta+num_theta*lphi])/w_norm << endl;
	};
      };
      // newline for every new angle \varphi
      file_spectrum_polar << endl;
    };
    // newline for every new energy
    file_spectrum_polar << endl;
    file_res << endl;
  }; // end of energy loop

  if (iv==1) cout << "Hasta la vista... " << endl;
  return 0;
};


