#ifndef TSURFF_SPECTRUM_HH
#define TSURFF_SPECTRUM_HH

#include <memory>
#include <vector>
#include <complex>
#include <grid.h>
#include <rwBinaryFile.hh>
#include <wavefunction.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <parameter.hh>
#ifdef HAVE_BOOST
#include <boost/timer.hpp>
#endif
#include <smallHelpers.hh>
#include <powers.hh>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

typedef std::complex<double> cplxd;
typedef std::unique_ptr<cplxd[]> cplxd_ptr;
typedef std::unique_ptr<double[]> double_ptr;

inline double Ylm_prefactor(long m) {
  if (m>=0)
    return 1.0;
  else if ((-m)%2==0)
    return 1.0;
  else
    return -1.0;
};

using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::ofstream;


/// This class saves the orbital and its derivatives at the t-SURFF boundary
///
class tsurffSaveWF {
  ofstream psi;
  ofstream dpsi_dr;
  long ind_R;
  // (let ((mylist '(("double" "delta-r" "param_ini") ("double" "R-tsurff" "param_tsurff") ("long" "qprop-dim" "param_ini") ("long" "ell-grid-size" "param_prop") ))) (param-list-for-class mylist))
  double delta_r;
  /// \f$\RI\f$
  double R_tsurff;
  /// \f$N_{\ell}\f$
  long ell_grid_size;
  long qprop_dim;
  grid the_grid;
public:
  /// Constructor
  /// @param_ini Parameters for the initial state
  /// @param_prop Parameters for the propagation
  tsurffSaveWF(const parameterListe &param_ini, const parameterListe &param_prop, const parameterListe &param_tsurff, grid prop_g, string filename="tsurff")
    : psi( filename+string("psi.raw"), std::ios::binary ), dpsi_dr( filename+string("-dpsidr.raw"), std::ios::binary ) {
    qprop_dim = param_ini.getLong("qprop-dim");
    delta_r = param_ini.getDouble("delta-r");
    R_tsurff = param_tsurff.getDouble("R-tsurff");
    ell_grid_size = param_prop.getLong("ell-grid-size");
    // the calculation of ind_R should be consistent with the rest of the program
    ind_R=prop_g.rindex(R_tsurff);
    // R_tsurff=prop_g.r(i_R);
    the_grid=prop_g;
  };
  /// Append the current wavefunction and its derivative to the binary files psi and dpsi\_dr
  ///
  /// qprop mode 34: \n
  /// \f$\text{psi} \leftarrow \phi_{\ell,m_0}(\RI) \quad \forall\ell\f$ \n
  /// \f$\text{dpsi\_dr} \leftarrow \partial_r\phi_{\ell,m_0}(r)|_{r=\RI} \quad \forall\ell\f$ \n
  /// qprop mode 44: \n
  /// \f$\text{psi} \leftarrow \phi_{\ell,m}(\RI) \quad \forall\ell,m\f$ \n
  /// \f$\text{dpsi\_dr} \leftarrow \partial_r\phi_{\ell,m}(r)|_{r=\RI} \quad \forall\ell,m\f$ \n
  void operator()(const wavefunction &wf) {
    const double one_over_12delta_r=1.0/(12.0*delta_r);
    const double two_over_3delta_r=2.0/(3.0*delta_r);
    if (qprop_dim==44) {
      for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	for (long i_m=-i_ell; i_m<i_ell+1; i_m++) {
	  const long ind_wave=the_grid.rlmindex(ind_R, i_ell, i_m);
	  const long ind_wave_n=the_grid.rlmindex(ind_R+1, i_ell, i_m);
	  const long ind_wave_p=the_grid.rlmindex(ind_R-1, i_ell, i_m);
	  const long ind_wave_n2=the_grid.rlmindex(ind_R+2, i_ell, i_m);
	  const long ind_wave_p2=the_grid.rlmindex(ind_R-2, i_ell, i_m);
	  psi.write(reinterpret_cast<const char*>(&wf[ind_wave]), sizeof(cplxd));
	  const cplxd deriv_wf = ((wf[ind_wave_n]-wf[ind_wave_p])*two_over_3delta_r-(wf[ind_wave_n2]-wf[ind_wave_p2])*one_over_12delta_r);
	  dpsi_dr.write(reinterpret_cast<const char*>(&deriv_wf), sizeof(cplxd));
	};
      };
    }
    else if (qprop_dim==34) {
      for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	const long ind_wave=the_grid.index(ind_R, i_ell, 0);
	const long ind_wave_n=the_grid.index(ind_R+1, i_ell, 0);
	const long ind_wave_p=the_grid.index(ind_R-1, i_ell, 0);
	const long ind_wave_n2=the_grid.index(ind_R+2, i_ell, 0);
	const long ind_wave_p2=the_grid.index(ind_R-2, i_ell, 0);
	psi.write(reinterpret_cast<const char*>(&wf[ind_wave]), sizeof(cplxd));
	const cplxd deriv_wf = ((wf[ind_wave_n]-wf[ind_wave_p])*two_over_3delta_r-(wf[ind_wave_n2]-wf[ind_wave_p2])*one_over_12delta_r);
	dpsi_dr.write(reinterpret_cast<const char*>(&deriv_wf), sizeof(cplxd));
      };
    }
    else {
      cerr << "No such option implemented for qprop dim " << qprop_dim << " in class tsurffSaveWF" << endl; 
    };
  };
};

/// Calculate PES for single active electron system interacting with a linearly polarized laser (\f$z\f$-direction) or a laser polarized in \f$xy\f$-plane. \n
/// Two methods for calculating the expansion of the probability amplitude for ionization are implemented. \n
/// In the case of expansion method 1 the coefficients of the expansion depend on the direction of the momentum
/// \f{equation*}{
/// b_1(\mathbf{k}) = \sum\limits^{N_{\ell}}_{\ell=0} \sum\limits^{\ell}_{m=-\ell} b_{1,\ell,m}(\mathbf{k}) Y_{\ell, m}(\Omega_k) \f}
/// wheras expansion method 2 is a "true" expansion in the sense that the expansion coefficients only depend on the absolute value of the momentum
/// \f{equation*}{
/// b_2(\mathbf{k}) = \sum\limits^{N_{\ell}}_{\ell=0} \sum\limits^{\ell}_{m=-\ell} b_{2,\ell,m}({k}) Y_{\ell, m}(\Omega_k) \f}
/// The template parameters vx, vy, vz are for classes which represent the components of the vector potential (See potentials.hh in one of the examples to learn how these should look like).
template<class vx, class vy, class vz>
class tsurffSpectrum {
  /// class which represents the \f$x\f$-component of the vector potential \f$A_x(t)\f$
  vx vecpot_x;
  /// class which represents the \f$y\f$-component of the vector potential \f$A_y(t)\f$
  vy vecpot_y;
  /// class which represents the \f$z\f$-component of the vector potential \f$A_z(t)\f$
  vz vecpot_z;

  ifstream dpsi_dr;
  ifstream psi;
  
  cplxd_ptr psi_surff;
  cplxd_ptr psi_A_surff;
  cplxd_ptr psi_A_cc_surff;
  cplxd_ptr psi_deriv_surff;
  cplxd_ptr psi_R;
  cplxd_ptr psi_R_deriv;
  double_ptr psi_buffer_R;
  double_ptr psi_buffer_R_deriv;
  cplxd_ptr partial_amplitude;
  // space for the momentum grid
  vector<double> thetas_surff;
  vector<double> phis_surff;
  vector<double> k_values;
  // (let ((mylist '(("double" "delta-r" "param_ini") ("double" "delta-t" "param_prop") ("double" "R-tsurff" "param_tsurff") ("double" "k-max-surff" "param_tsurff") ("long" "num-k-surff" "param_tsurff") ("long" "num-theta-surff" "param_tsurff") ("long" "num-phi-surff" "param_tsurff") ("long" "delta-k-scheme" "param_tsurff") ("long" "cache-size-t" "param_tsurff") ("long" "expansion-scheme" "param_tsurff") ("long" "qprop-dim" "param_ini") ("long" "ell-grid-size" "param_prop") ("long" "initial-m" "param_ini")))) (param-list-for-class mylist))

  long expansion_scheme;
  double delta_r;
  double delta_t;
  /// \f$\RI\f$
  double R_tsurff;
  double k_max_surff;
  long num_k_surff;
  long num_theta_surff;
  long num_phi_surff;
  long delta_k_scheme;
  long cache_size_t;
  long qprop_dim;
  /// \f$N_{\ell}\f$
  long ell_grid_size;
  long initial_m;

  long ell_m_grid_size;
  
  double duration;

  int i_proc;
  int num_proc;
  long  num_k_proc;
public:
  tsurffSpectrum(const parameterListe &param_ini, const parameterListe &param_prop, const parameterListe &param_tsurff, vx vecpotx, vy vecpoty, vz vecpotz, string filename="tsurff")
    : thetas_surff(101, 0.0), phis_surff(101, 0.0), k_values(100, 0.0), psi(filename+string("psi.raw"), std::ios::binary), dpsi_dr(filename+string("-dpsidr.raw"), std::ios::binary), vecpot_x(vecpotx), vecpot_y(vecpoty), vecpot_z(vecpotz) {
    expansion_scheme = param_tsurff.getLong("expansion-scheme");
    delta_r = param_ini.getDouble("delta-r");
    delta_t = param_prop.getDouble("delta-t");
    R_tsurff = param_tsurff.getDouble("R-tsurff");
    k_max_surff = param_tsurff.getDouble("k-max-surff");
    num_k_surff = param_tsurff.getLong("num-k-surff");
    num_theta_surff = param_tsurff.getLong("num-theta-surff");
    num_phi_surff = param_tsurff.getLong("num-phi-surff");
    delta_k_scheme = param_tsurff.getLong("delta-k-scheme");
    cache_size_t = param_tsurff.getLong("cache-size-t");
    qprop_dim = param_ini.getLong("qprop-dim");
    ell_grid_size = param_prop.getLong("ell-grid-size");
    initial_m = param_ini.getLong("initial-m");

    // size of grid for ells and ms: ell^2 for dim 44 and ell for 34
    ell_m_grid_size=(qprop_dim==34)?ell_grid_size:ell_grid_size*ell_grid_size;
    
    // the calculation of ind_R should be consistent with the rest of the program
    grid prop_g(qprop_dim, 0, delta_r, 0.0); // this grid instance is only used to calculate an r-index and r from index
    const long ind_R=prop_g.rindex(R_tsurff);
    // Caution: value from parameter file is overwritten by value on grid point
    R_tsurff=prop_g.r(ind_R);

    // This might be a bug waiting to happen!!
    const double pulse_duration=(qprop_dim==34)?vecpot_z.get_duration():vecpot_x.get_duration();
    
    // how long do the slowest electrons have time to reach the t-SURFF boundary
    const double time_surff=param_tsurff.getDouble("R-tsurff")/param_tsurff.getDouble("p-min-tsurff");
    duration=pulse_duration+time_surff;

    i_proc=0;
    num_proc=1;
    num_k_proc=num_k_surff;
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &i_proc);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    num_k_proc=num_k_surff/num_proc;
    // last proc has to work more sometimes..
    if (i_proc==(num_proc-1))
      num_k_proc+=num_k_surff%num_proc;
    cout << " I am process " << i_proc << " of " << num_proc << " processes. My share of k values: " << num_k_proc << endl;
#endif
    
    prep_mom_grid();
    prep_space();
  };

  /// choose the correct function for the necessary time integration depending on polarization and expansion method \n
  /// read wavefunction data from disk \n
  /// expansion_method=1, qprop_dim=34: \n
  /// call time_integration_theta_phi(long, long) \n
  /// expansion_method=1, qprop_dim=44: \n
  /// call time_integration_theta_phi_ell_m(long, long) \n
  /// expansion_method=2, qprop_dim=34: \n
  /// call time_integration_ell1_ell2(long, long) \n
  /// expansion_method=2, qprop_dim=44: \n
  /// call time_integration_ell1_m1_ell2_m2(long, long) \n
  void time_integration() {
    long entries_read(0);
    long all_entries_read(0);
    long offset(0);
    long start_t(0);
    do {
#ifdef HAVE_BOOST
      boost::timer tim;
#endif
      // ell_m_grid_size * cache_size_t entries are read
      entries_read=read_psi_data()/ell_m_grid_size;
      if (expansion_scheme==1 && qprop_dim==34)
	time_integration_theta_phi(start_t, start_t+entries_read);
      else if (expansion_scheme==2 && qprop_dim==34) {
	time_integration_ell1_ell2(start_t, start_t+entries_read);
      }
      else if (expansion_scheme==2 && qprop_dim==44) {
	time_integration_ell1_m1_ell2_m2(start_t, start_t+entries_read);
      }
      else if (expansion_scheme==1 && qprop_dim==44) {
	time_integration_theta_phi_ell_m(start_t, start_t+entries_read);
      }
      else {
	cerr << "combination of expansion scheme " << expansion_scheme << " and qprop dimension " << qprop_dim << " is not implemented." << endl;
	exit(-1);
      };
      start_t+=entries_read;
      all_entries_read+=entries_read;
      if (i_proc==0)
	cout << "timestep " << all_entries_read << " of " << long(duration/delta_t)+1 << endl;
      offset+=cache_size_t;
#ifdef HAVE_BOOST
      cout << "time for call of time_integration() per time step:  " << tim.elapsed()/double(entries_read) << " for proc " << i_proc << endl;
#endif
    } while (entries_read==cache_size_t);
  };
  
  /// choose the correct function for calculating the PES depending on polarization and expansion method
  /// \par expansion_method=1, qprop_dim=34: \n
  /// call eval_spectrum_34_theta_phi() \n
  /// \par expansion_method=1, qprop_dim=44: \n
  /// call eval_spectrum_34_theta_phi() \n
  /// \par expansion_method=2, qprop_dim=34: \n
  /// call eval_partial_spectra_34_ell1_ell2() \n
  /// call eval_spectrum_34_ell1_ell2() \n
  /// \par expansion_method=2, qprop_dim=44: \n
  /// call eval_partial_spectra_44_ell1_ell2() \n
  /// call eval_spectrum_44_ell1_ell2() \n
  void polar_spectrum() {
    if (expansion_scheme==1 && qprop_dim==34) {
      eval_spectrum_34_theta_phi();
    }
    else if (expansion_scheme==1 && qprop_dim==44) {
      eval_spectrum_44_theta_phi();
    }
    else if (expansion_scheme==2 && qprop_dim==34) {
      eval_partial_spectra_34_ell1_ell2();
      eval_spectrum_34_ell1_ell2();
    }
    else if (expansion_scheme==2 && qprop_dim==44) {
      eval_partial_spectra_44_ell1_ell2();
      eval_spectrum_44_ell1_ell2();
    }
    else {
      cerr << "expansion scheme " << expansion_scheme << " not implemented" << endl;
      exit(-1);
    };
  };

  /// print the partial amplitudes to a file if expansion method is 2
  /// and do nothing otherwise
  void print_partial_amplitudes() {
    if (expansion_scheme==2) {
      ofstream debug_dat("tsurff-partial"+to_string(i_proc)+".dat");
      if (qprop_dim==34 && expansion_scheme==2) {
	for (long i_k=0; i_k<num_k_proc; i_k++) {
	  double total_spectrum_k(0.0);
	  const double k=k_values[i_k];
	  debug_dat << k*k*0.5 << " " << k;
	  for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	    const long ind_ampl=index_ell_k(i_ell, i_k);
	    // debug_dat << " " << real(partial_amplitude[index_ell_k(i_ell, i_k)]) << " " << imag(partial_amplitude[index_ell_k(i_ell, i_k)]);
	    const double norm_a_ell_k=norm(partial_amplitude[ind_ampl]);
	    total_spectrum_k+=norm_a_ell_k;
	    debug_dat << " " << norm_a_ell_k*k;
	  };
	  debug_dat << " " << total_spectrum_k*k << endl;
	};
      }
      else if (qprop_dim==44 && expansion_scheme==2) {
	for (long i_k=0; i_k<num_k_proc; i_k++) {
	  double total_spectrum_k(0.0);
	  const double k=k_values[i_k];
	  debug_dat << k*k*0.5 << " " << k;
	  for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	    for (long i_m=-(i_ell); i_m<(i_ell+1); i_m++) {
	      const long ind_ampl=index_ell_m_k(i_ell, i_m, i_k);
	      // debug_dat << " " << real(partial_amplitude[ind_ampl]) << " " << imag(partial_amplitude[ind_ampl]);
	      const double norm_a_ell_m_k=norm(partial_amplitude[ind_ampl]);
	      total_spectrum_k+=norm_a_ell_m_k;
	      debug_dat << " " << norm_a_ell_m_k*k;
	    };
	  };
	  debug_dat << " " << total_spectrum_k*k << endl;
	};
      };
    };
  };
  
private:
  /// setup the grid in momentum space using the parameters from tsurff.param
  /// side effects: writes to k_values, thetas_surff, phis_surff
  void prep_mom_grid() {
    if (num_k_proc>100)
      k_values.resize(num_k_proc);
    long i_k_min=0;
    long i_k_max=num_k_surff;
#ifdef HAVE_MPI
    i_k_min=i_proc * ( num_k_surff / num_proc );
    i_k_max=(i_proc+1) * ( num_k_surff / num_proc );
    if (i_proc==(num_proc-1))
      i_k_max+=num_k_surff%num_proc;
#endif
    if (delta_k_scheme==1) {
      const double delta_k_surff=k_max_surff/double(num_k_surff);
      for (long i_loop=i_k_min; i_loop<i_k_max; i_loop++) {
	k_values[i_loop-i_k_min]=double(i_loop)*delta_k_surff;
      };
    }
    else if (delta_k_scheme==2) {
      const double delta_k_surff=pow2(k_max_surff)/double(num_k_surff);
      for (long i_loop=i_k_min; i_loop<i_k_max; i_loop++) {
	k_values[i_loop-i_k_min]=sqrt(double(i_loop)*delta_k_surff);
      };
    }
    else {
      cerr << "delta_k_scheme can only be 1 for equal spaceing in k or 2 for equal spacing in k^2" << endl;
      exit(-42);
    };
    cout << " I am process " << i_proc << " My first k value " << k_values[0] << endl;
    // min and max \theta should always be 0 and \pi
    // it is more convenient if number of thetas is odd
    if (num_theta_surff%2==0) num_theta_surff++;
    if (num_theta_surff<3) num_theta_surff=3;
    if (num_theta_surff>101) thetas_surff.resize(num_theta_surff);
    const double delta_theta=M_PI/double(num_theta_surff-1);

    for (long i_loop=0; i_loop<num_theta_surff; i_loop++) {
      thetas_surff[i_loop]=(i_loop)*delta_theta;
    };
    
    // min and max \phi should always be 0 and 2 \pi
    // num_phi_surff=50;
    // it is more convenient if number of phis is even
    // if (num_phi_surff%2==1) num_theta_surff++;
    // if (num_phi_surff<4) num_phi_surff=4;
    if (num_phi_surff>101) phis_surff.resize(num_phi_surff);
    const double delta_phi=2.0*M_PI/double(num_phi_surff);

    for (long i_loop=0; i_loop<num_phi_surff; i_loop++) {
      phis_surff[i_loop]=(i_loop)*delta_phi;
    };
  };

  /// allocate space for the time integrals and PES results
  void prep_space() {
    // the size of the result after summing over the $\ell$s and $m$s
    // const long num_spectrum_size_surff=num_k_proc*num_theta_surff*num_phi_surff;
    // spectrum_surff.reset(new cplxd[num_spectrum_size_surff]);
   
    // now we need two chunks of memory for the time integral of the orbitals (space for every theta, phi, k, l, m for scheme 1 or k, l1, l2, m1, m2 for scheme 2)
    // more space is needed for the 44 mode
    const long num_psi_size_surff=(expansion_scheme==1)?(ell_m_grid_size*num_theta_surff*num_phi_surff*num_k_proc):ell_m_grid_size*num_k_proc*ell_m_grid_size;

    if (expansion_scheme==2)
      partial_amplitude.reset(new cplxd[ell_m_grid_size*num_k_proc]);
    
    psi_surff.reset(new cplxd[num_psi_size_surff]);
    psi_A_surff.reset(new cplxd[num_psi_size_surff]);
    psi_deriv_surff.reset(new cplxd[num_psi_size_surff]);
    // To be on the save side initialize with zeros
    for (long i_loop=0; i_loop<num_psi_size_surff; i_loop++) {
      psi_surff[i_loop]=cplxd(0.0);
      psi_A_surff[i_loop]=cplxd(0.0);
      psi_deriv_surff[i_loop]=cplxd(0.0);
    };
    // this is only used in 44 mode
    if (qprop_dim==44) {
      psi_A_cc_surff.reset(new cplxd[num_psi_size_surff]);
      for (long i_loop=0; i_loop<num_psi_size_surff; i_loop++) {	
	psi_A_cc_surff[i_loop]=cplxd(0.0);
      };
    };
    // space for loading wf data
    psi_R.reset(new cplxd[ell_m_grid_size*cache_size_t]);
    psi_R_deriv.reset(new cplxd[ell_m_grid_size*cache_size_t]);
    psi_buffer_R.reset(new double[ell_m_grid_size*cache_size_t*2]);
    psi_buffer_R_deriv.reset(new double[ell_m_grid_size*cache_size_t*2]);
    cout << "saving the time integrals will take " << double(((qprop_dim==44)?4:3)*num_psi_size_surff*16)/1.e6 << "MB of space" << endl;
  };

  // First process reads the data from file and broadcasts it to the rest. The buffers are necessary because not all mpi implementations know copmlex<double>.
  // The chunk of data read in one step is determined by cache_size_t (a number of timesteps).
  long read_psi_data() {
#ifdef HAVE_BOOST
    boost::timer tim;
#endif
    // Actually MPI_Bcast will convert this back to int..what a waste
    long count=0;
    if (i_proc==0) {
      while (psi.good() && count<(cache_size_t*ell_m_grid_size)) { 
	psi.read(reinterpret_cast< char* >(&psi_R[count]), sizeof(cplxd));
	dpsi_dr.read(reinterpret_cast< char* >(&psi_R_deriv[count]), sizeof(cplxd));
#ifdef HAVE_MPI
	psi_buffer_R[2*count]=real(psi_R[count]); psi_buffer_R[2*count+1]=imag(psi_R[count]);
	psi_buffer_R_deriv[2*count]=real(psi_R_deriv[count]); psi_buffer_R_deriv[2*count+1]=imag(psi_R_deriv[count]);
#endif
	// check if there is more data or end of file
	psi.peek();
	count++;      
      };
    };
#ifdef HAVE_MPI
    MPI_Bcast(&count, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(psi_buffer_R_deriv.get(), 2*count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(psi_buffer_R.get(), 2*count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (long i=0; i<count; i++) {
      psi_R[i]=cplxd(psi_buffer_R[2*i], psi_buffer_R[2*i+1]);
      psi_R_deriv[i]=cplxd(psi_buffer_R_deriv[2*i], psi_buffer_R_deriv[2*i+1]);
    };
#endif

#ifdef HAVE_BOOST
    cout << "time for data reading: " << tim.elapsed() << " for proc "<< i_proc << endl;
#endif

    return count;
  };

  /// calculate an index (useful in expansion mode 1 and qprop dim 34)
  ///
  /// \f$i=i_{\ell} + i_{\theta} N_{\ell} + i_{k} N_{\ell} N_{\theta}\f$ 
  inline long index_ell_theta_k(long i_ell, long i_theta, long i_k) const {
    return  i_ell + i_theta*ell_grid_size + i_k*ell_grid_size*num_theta_surff;
  };

  /// Perform the necessary time integrations for linear polarization and expansion method 1
  ///
  /// side-effects: writes to psi_surff, psi_A_surff, psi_deriv_surff
  /// \f{align*}{
  /// \text{psi\_A\_surff} \leftarrow I_{2,\ell,m}(k,k_z,\RI,T) &=\int_0^T dt H(t) e^{i k^2 t/2 + ik_z \alpha_z(t)} A_z(t) \phi_{\ell, m}(\RI,t) \\
  /// \text{psi\_surff} \leftarrow I_{3,\ell,m}(k,k_z,\RI,T) &=\int_0^T dt H(t) e^{i k^2 t/2 + ik_z \alpha_z(t)} \phi_{\ell, m}(\RI,t) \\
  /// \text{psi\_deriv\_surff} \leftarrow I_{4,\ell,m}(k,k_z,\RI,T) &=\int_0^T dt H(t) e^{i k^2 t/2 + ik_z \alpha_z(t)} \partial_r \phi_{\ell, m}(r,t) \vert_{r=\RI} \\
  /// \f}
  /// The Hanning window \f$H(t)\f$ is implemented in the function hanning_win(double) const.
  /// The vector potential \f$A_z(t)\f$ and the time integral of the vector potential \f$\alpha_z(t)\f$
  /// are implemented in member functions of the tsurffSpectrum::vecpot_z object
  void time_integration_theta_phi(long start_t, long stop_t) {
    // time integration for t-surff
    for (long i_k=0; i_k<num_k_proc; i_k++) {
      const double k=k_values[i_k];
      for (long i_theta=0; i_theta<num_theta_surff; i_theta++) {
	const double cos_theta=cos(thetas_surff[i_theta]);
	for (long i_t=start_t; i_t<stop_t; i_t++) {
	  const long cache_i_t=i_t-start_t;
	  const double time=delta_t*double(i_t);
	  const double hanning_t=hanning_win(time);
	  const double alpha_t=vecpot_z.integral(time);
	  const double A_t=vecpot_z(time, 0);
	  // H(t) e^{i k^2 t/2 + ik_z \alpha_z(t)}
	  const cplxd hanning_t_exp_volkov=hanning_t*exp(cplxd(0., 1.)*k*k*time*0.5 + cplxd(0., 1.)*k*cos_theta*alpha_t);
	  for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	    const long ind_surff=index_ell_theta_k(i_ell, i_theta, i_k);
	    const long ind_wave=i_ell + cache_i_t*ell_grid_size;
	    psi_surff[ind_surff]      +=hanning_t_exp_volkov*psi_R[ind_wave];
	    psi_A_surff[ind_surff]    +=A_t*hanning_t_exp_volkov*psi_R[ind_wave];
	    psi_deriv_surff[ind_surff]+=hanning_t_exp_volkov*psi_R_deriv[ind_wave];
	  };
	};
      };
    };
  };

  /// calculate an index (useful in expansion mode 1 and qprop dim 44)
  ///
  /// \f$i=(i_{\ell} + 1)^2 - (i_{\ell} + 1) + i_{\theta} N_{\ell}^2 + i_{\phi} N_{\ell}^2 N_{\theta} + i_{k} N_{\ell}^2 N_{\theta} N_{\phi} \f$ 
  inline long index_ell_m_theta_phi_k(long i_ell, long i_m, long i_theta, long i_phi, long i_k) const {
    return  ((i_ell+1)*(i_ell+1)-(i_ell+1)+i_m) + i_theta*ell_m_grid_size + i_phi*num_theta_surff*ell_m_grid_size + i_k*num_phi_surff*num_theta_surff*ell_m_grid_size;
  };

  /// Perform the necessary time integrations for polarization in the \f$xy\f$-plane and expansion method 1
  ///
  /// side-effects: writes to psi_surff, psi_A_surff, psi_A_cc_surff,psi_deriv_surff
  /// \f{align*}{
  /// \text{psi\_A\_cc\_surff} \leftarrow I_{0,\ell,m}(k,k_x,k_y,\RI,T) &=\int_0^T dt H(t) e^{i k^2 t/2 + i k_x \alpha_x(t) + i k_y \alpha_y(t)} \tilde{A}^{*}(t) \phi_{\ell, m}(\RI,t) \\
  /// \text{psi\_A\_surff} \leftarrow I_{1,\ell,m}(k,k_x,k_y,\RI,T) &=\int_0^T dt H(t) e^{i k^2 t/2 + i k_x \alpha_x(t) + i k_y \alpha_y(t)} \tilde{A}(t) \phi_{\ell, m}(\RI,t) \\
  /// \text{psi\_surff} \leftarrow I_{3,\ell,m}(k,k_x,k_y,\RI,T) &=\int_0^T dt H(t) e^{i k^2 t/2 + i k_x \alpha_x(t) + i k_y \alpha_y(t)} \phi_{\ell, m}(\RI,t) \\
  /// \text{psi\_deriv\_surff} \leftarrow I_{4,\ell,m}(k,k_x,k_y,\RI,T) &=\int_0^T dt H(t) e^{i k^2 t/2 + i k_x \alpha_x(t) + i k_y \alpha_y(t)} \partial_r \phi_{\ell, m}(r,t) \vert_{r=\RI} \\
  /// \f}
  /// where \f$\tilde{A}=A_x+i A_y\f$ \n
  /// The Hanning window \f$H(t)\f$ is implemented in the function hanning_win(double) const.
  /// The vector potential components \f$A_x(t)\f$, \f$A_y(t)\f$ and the time integrals of the vector potential \f$\alpha_x(t)\f$, \f$\alpha_y(t)\f$
  /// are implemented in member functions of the tsurffSpectrum::vecpot_x and tsurffSpectrum::vecpot_y objects.
  void time_integration_theta_phi_ell_m(long start_t, long stop_t) {
    // time integration for t-surff
    for (long i_k=0; i_k<num_k_proc; i_k++) {
      const double k=k_values[i_k];
      for (long i_phi=0; i_phi<num_phi_surff; i_phi++) {
	const double cos_phi=cos(phis_surff[i_phi]);
	const double sin_phi=sin(phis_surff[i_phi]);
	for (long i_theta=0; i_theta<num_theta_surff; i_theta++) {
	  const double sin_theta=sin(thetas_surff[i_theta]);
	  for (long i_t=start_t; i_t<stop_t; i_t++) {
	    const long cache_i_t=i_t-start_t;
	    const double time=delta_t*double(i_t);
	    const double hanning_t=hanning_win(time);
	    const double alpha_x_t=vecpot_x.integral(time);
	    const double alpha_y_t=vecpot_y.integral(time);
	    const cplxd A_tilde(vecpot_x(time, 0), vecpot_y(time, 0));
	    const cplxd hanning_t_exp_volkov=hanning_t*exp(cplxd(0., 1.)*k*k*time*0.5 + cplxd(0., 1.)*k*cos_phi*sin_theta*alpha_x_t + cplxd(0., 1.)*sin_theta*k*sin_phi*alpha_y_t);
	    for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	      for (long i_m=-i_ell; i_m<(i_ell+1); i_m++) {
		const long ind_surff=index_ell_m_theta_phi_k(i_ell, i_m, i_theta, i_phi, i_k);
		const long ind_wave=index_ell1_m1_i_t(i_ell, i_m, cache_i_t);
		psi_surff[ind_surff]      +=hanning_t_exp_volkov*psi_R[ind_wave];
		psi_A_surff[ind_surff]    +=A_tilde*hanning_t_exp_volkov*psi_R[ind_wave];
		psi_A_cc_surff[ind_surff] +=conj(A_tilde)*hanning_t_exp_volkov*psi_R[ind_wave];
		psi_deriv_surff[ind_surff]+=hanning_t_exp_volkov*psi_R_deriv[ind_wave];
	      };
	    };
	  };
	};
      };
    };
  };

  /// calculate an index (useful in expansion mode 2 and qprop dim 34)
  ///
  /// \f$i=i_{\ell_1} + i_{\ell_2} N_{\ell} + i_{k} N_{\ell}^2 \f$ 
  inline long index_ell1_ell2_k(long i_ell1, long i_ell2, long i_k) const {
    return  i_ell1 + i_ell2*ell_grid_size + i_k*ell_grid_size*ell_grid_size;
  };

  /// Perform the necessary time integrations for linear polarization and expansion method 2
  ///
  /// side-effects: writes to psi_surff, psi_A_surff, psi_deriv_surff
  /// \f{align*}{
  /// \text{psi\_A\_surff} \leftarrow I_{2,\ell_1,m_1,\ell_2}(k) &= \int_0^T \!dt H(t)   e^{ik^2t/2} j_{\ell_2}(k \alpha(t)) A_z(t) \phi_{\ell_1,m_1} \\
  /// \text{psi\_surff} \leftarrow I_{3,\ell_1,m_1,\ell_2}(k) &= \int_0^T \!dt H(t)   e^{ik^2t/2} j_{\ell_2}(k \alpha(t)) \phi_{\ell_1,m_1} \\ 
  /// \text{psi\_deriv\_surff} \leftarrow I_{4,\ell_1,m_1,\ell_2}(k) &= \int_0^T \!dt H(t)   e^{ik^2t/2} j_{\ell_2}(k \alpha(t)) \partial_r \phi_{\ell_1,m_1} \\
  /// \f}
  /// The Hanning window \f$H(t)\f$ is implemented in the function hanning_win(double) const.
  /// The vector potential \f$A_z(t)\f$ and the time integral of the vector potential \f$\alpha_z(t)\f$
  /// are implemented in member functions of the tsurffSpectrum::vecpot_z object.
  /// The function
  /// fill_bessel_array(double, double, double*) handles the computation of the Bessel functions \f$j_{\ell_2}(k\alpha(t))\f$.
  void time_integration_ell1_ell2(long start_t, long stop_t) {
    // time integration for t-surff
    double j_l2kalpha[ell_grid_size+1];
    for (long i_k=0; i_k<num_k_proc; i_k++) {
      const double k=k_values[i_k];
	for (long i_t=start_t; i_t<stop_t; i_t++) {
	  const long cache_i_t=i_t-start_t;
	  const double time=delta_t*double(i_t);
	  const double hanning_t=hanning_win(time);
	  const double alpha_t=vecpot_z.integral(time);
	  const double A_t=vecpot_z(time, 0);
	  // e^{ik^2t/2} j_{\ell_2}(k \alpha(t))
	  const cplxd hanning_t_exp = hanning_t*exp(cplxd(0., 1.)*k*k*time*0.5);
	  fill_bessel_array(k, alpha_t, j_l2kalpha);
	  for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
	    for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	      const long ind_surff=index_ell1_ell2_k(i_ell, i_ell2, i_k);
	      const long ind_wave=i_ell + cache_i_t*ell_grid_size;
	      psi_surff[ind_surff] += hanning_t_exp*j_l2kalpha[i_ell2]*psi_R[ind_wave];
	      psi_A_surff[ind_surff] += A_t*hanning_t_exp*j_l2kalpha[i_ell2]*psi_R[ind_wave];
	      psi_deriv_surff[ind_surff] += hanning_t_exp*j_l2kalpha[i_ell2]*psi_R_deriv[ind_wave];
	  };
	};
      };
    };
  };

  inline long index_ell1_m1_ell2_m2_k(long i_ell1, long i_m1, long i_ell2, long i_m2, long i_k) const {
    const long ell_grid_size_2=ell_grid_size*ell_grid_size;
    return  ((i_ell1+1)*(i_ell1+1)-(i_ell1+1)+i_m1) + ((i_ell2+1)*(i_ell2+1)-(i_ell2+1)+i_m2)*ell_grid_size_2 + i_k*ell_grid_size_2*ell_grid_size_2;
  };

  inline long index_ell1_m1_i_t(long i_ell1, long i_m1, long i_t) const {
    const long ell_grid_size_2=ell_grid_size*ell_grid_size;
    return  ((i_ell1+1)*(i_ell1+1)-(i_ell1+1)+i_m1) + i_t*ell_grid_size_2;
  };

  /// Perform the necessary time integrations for polarization in \f$xy\f$-plane and expansion method 2
  ///
  /// side-effects: writes to psi_surff, psi_A_surff, psi_A_cc_surff, psi_deriv_surff
  /// \f{align*}{
  /// \text{psi\_A\_cc\_surff} \leftarrow I_{0,\ell_1,m_1,\ell_2,m_2}(k) &= \int_0^T \!dt H(t)  e^{ik^2t/2} j_{\ell_2}(k \alpha(t)) {(-1)}^{m_2} Y^*_{\ell_2,-m_2}(\Omega_{\alpha}) \tilde{A}^*(t) \phi_{\ell_1,m_1} \\
  ///  \text{psi\_A\_surff} \leftarrow I_{1,\ell_1,m_1,\ell_2,m_2}(k) &= \int_0^T \!dt H(t) e^{ik^2t/2} j_{\ell_2}(k \alpha(t)) {(-1)}^{m_2} Y^*_{\ell_2,-m_2}(\Omega_{\alpha}) \tilde{A}(t) \phi_{\ell_1,m_1} \\
  /// \text{psi\_surff} \leftarrow I_{3,\ell_1,m_1,\ell_2,m_2}(k) &= \int_0^T \!dt H(t) e^{ik^2t/2} j_{\ell_2}(k \alpha(t)) {(-1)}^{m_2} Y^*_{\ell_2,-m_2}(\Omega_{\alpha}) \phi_{\ell_1,m_1} \\
  /// \text{psi\_deriv\_surff} \leftarrow I_{4,\ell_1,m_1,\ell_2,m_2}(k) &= \int_0^T \!dt H(t) e^{ik^2t/2} j_{\ell_2}(k \alpha(t)) {(-1)}^{m_2} Y^*_{\ell_2,-m_2}(\Omega_{\alpha}) \partial_r \phi_{\ell_1,m_1}
  /// \f}
  /// where \f$\tilde{A}=A_x+i A_y\f$ \n
  /// The Hanning window \f$H(t)\f$ is implemented in the function hanning_win(double) const.
  /// The vector potential components \f$A_x(t)\f$, \f$A_y(t)\f$ and the time integrals of the vector potential \f$\alpha_x(t)\f$, \f$\alpha_y(t)\f$
  /// are implemented in member functions of the tsurffSpectrum::vecpot_x and tsurffSpectrum::vecpot_y objects.
  /// The function
  /// fill_bessel_array(double, double, double*) handles the computation of the Bessel functions \f$j_{\ell_2}(k\alpha(t))\f$.
  void time_integration_ell1_m1_ell2_m2(long start_t, long stop_t) {
    // time integration for t-surff
    double j_l2kalpha[ell_grid_size+1];
    for (long i_k=0; i_k<num_k_proc; i_k++) {
      const double k=k_values[i_k];
      for (long i_t=start_t; i_t<stop_t; i_t++) {
	const long cache_i_t=i_t-start_t;
	const double time=delta_t*double(i_t);
	const double hanning_t=hanning_win(time);
	const double alpha_x=vecpot_x.integral(time);
	const double alpha_y=vecpot_y.integral(time);
	const double abs_alpha=sqrt(pow2(alpha_x)+pow2(alpha_y));
	const double phi_alpha=atan2(alpha_y, alpha_x);
	const cplxd A_tilde(vecpot_x(time, 0), vecpot_y(time, 0));
	// const cplxd A_tilde(0,0);
	// e^{ik^2t/2} j_{\ell_2}(k \alpha(t))
	const cplxd hanning_t_exp = hanning_t*exp(cplxd(0., 1.)*k*k*time*0.5);
	fill_bessel_array(k, abs_alpha, j_l2kalpha);
	for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
	  for (long i_m2=-i_ell2; i_m2<i_ell2+1; i_m2++) {
	    const cplxd Y_ell2_m2=Ylm_prefactor(i_m2)*gsl_sf_legendre_sphPlm(i_ell2, labs(i_m2), 0)*exp(cplxd(0, double(i_m2)*phi_alpha));
	    const cplxd A_tilde_hanning_t_exp_j_l2kalpha_Y_ell2_m2=A_tilde*hanning_t_exp*j_l2kalpha[i_ell2]*Y_ell2_m2;
	    const cplxd A_tilde_cc_hanning_t_exp_j_l2kalpha_Y_ell2_m2=conj(A_tilde)*hanning_t_exp*j_l2kalpha[i_ell2]*Y_ell2_m2;
	    const cplxd hanning_t_exp_j_l2kalpha_Y_ell2_m2=hanning_t_exp*j_l2kalpha[i_ell2]*Y_ell2_m2;	    
	    for (long i_ell1=0; i_ell1<ell_grid_size; i_ell1++) {
	      for (long i_m1=-i_ell1; i_m1<i_ell1+1; i_m1++) {
		const long ind_surff=index_ell1_m1_ell2_m2_k(i_ell1, i_m1, i_ell2, i_m2, i_k);
		const long ind_wave=index_ell1_m1_i_t(i_ell1, i_m1, cache_i_t);
		psi_surff[ind_surff] += hanning_t_exp_j_l2kalpha_Y_ell2_m2*psi_R[ind_wave];
		psi_A_surff[ind_surff] += A_tilde_hanning_t_exp_j_l2kalpha_Y_ell2_m2*psi_R[ind_wave];
		psi_A_cc_surff[ind_surff] += A_tilde_cc_hanning_t_exp_j_l2kalpha_Y_ell2_m2*psi_R[ind_wave];
		psi_deriv_surff[ind_surff] += hanning_t_exp_j_l2kalpha_Y_ell2_m2*psi_R_deriv[ind_wave];
		
		// psi_surff[ind_surff] += hanning_t_exp*j_l2kalpha[i_ell2]*Y_ell2_m2*psi_R[ind_wave];
		// psi_A_surff[ind_surff] += A_tilde*hanning_t_exp*j_l2kalpha[i_ell2]*Y_ell2_m2*psi_R[ind_wave];
		// psi_A_cc_surff[ind_surff] += conj(A_tilde)*hanning_t_exp*j_l2kalpha[i_ell2]*Y_ell2_m2*psi_R[ind_wave];
		// psi_deriv_surff[ind_surff] += hanning_t_exp*j_l2kalpha[i_ell2]*Y_ell2_m2*psi_R_deriv[ind_wave];
	      };
	    };
	  };
	};
      };
    };
  };

  /// half of a hanning window
  ///
  /// \f{equation*}{
  /// H(t)=
  ///   \begin{cases}
  ///    1 & \text{if } t < T/2 \\
  ///    [1-\cos(2\pi t/T)]/2       & \text{if } t \geq T/2 .
  ///   \end{cases}
  /// \f}
  inline double hanning_win(double time) const {
    double hanning = 1.0;
    if (time>(0.5*(duration)))
      hanning = 0.5*(1.0-cos(2.0*M_PI*time/(duration)));
    return hanning;
  };

  /// calculate Bessel functions for the argument \f$kR\f$ up to \f$N_{\ell}\f$
  inline void fill_bessel_array(double k, double R, double* j_lkr) {
    // Function: double gsl_sf_bessel_jl (int l, double x)
    // Function: int gsl_sf_bessel_jl_e (int l, double x, gsl_sf_result * result)
    // These routines compute the regular spherical Bessel function of order l, j_l(x), for l>= 0 and x >= 0.
    // calculate Bessel functions for increasing \ell until |j_{\ell}(kR)|<eps_j
    const double eps_j(1.0e-100);
    long max_ell_bessel(ell_grid_size+1);
    const double a_sign=(k*R<0.0)?-1.0:1.0;
    long i_ell1(0);
    do {
      j_lkr[i_ell1]=gsl_sf_bessel_jl(i_ell1, a_sign*k*R);
      i_ell1++;
    } while (fabs(j_lkr[i_ell1-1])>eps_j && i_ell1<max_ell_bessel);
    if (i_ell1==0) {
      exit(-1);
      cerr << "Bessel calculation has failed!" << endl;
    };
    // now set the rest of the ells to zero
    for (long i_ell=i_ell1; i_ell<ell_grid_size+1; i_ell++) {
      j_lkr[i_ell]=0.0;
    };
    // treat negative arguments: sign(k*R)^{\ell} j(|k R|)_{\ell}
    if (a_sign<0.0) {
      for (long i_ell=1; i_ell<ell_grid_size+1; i_ell+=2) {
	j_lkr[i_ell]*=-1.0;
      };
    };      
  };
  
  /// calculate angle and energy resolved spectrum for linear polarization and expansion method 1 \n
  /// side effects: writes to file tsurff-polari_proc.dat \n
  /// reads time integrals from arrays psi_A_surff psi_surff and psi_deriv_surff \n
  /// dependencies: Ylm_prefactor and gsl_sf_legendre_sphPlm for calculating the spherical harmonics \f$Y_{\ell m}\f$,
  ///               fill_bessel_array for calculating spherical bessel functions \n
  /// The probability amplitude is
  /// \f{equation*}{
  ///   \label{eq:pol-tsurff-ell-m-sum}
  /// b_1(\mathbf{k}) = \sum\limits^{N_{\ell}}_{\ell=0} b_{1,\ell,m}(\mathbf{k}) Y_{\ell, m}(\Omega_k)
  /// \f}
  /// where
  ///   \f{multline*}{
  ///   b_{1,\ell, m}(\mathbf{k})
  ///   = \frac{1}{(2\pi)^{1/2}} \RI
  ///   \Bigg[ 2  c_{\ell-1, m} (-i)^{\ell} j_{\ell}(k \RI) I_{2,\ell-1} 
  ///   +  2  c_{\ell, m} (-i)^{\ell} j_{\ell}(k \RI)  I_{2,\ell+1} \\
  ///   +  (-i)^{\ell+1} j_{\ell}(k \RI)  
  ///   \left( -\RI^{-1} I_{3,\ell, m} +  I_{4,\ell, m} \right)
  ///   + (-i)^{\ell+1} k\left(j_{\ell+1}(k\RI) - \frac{\ell}{k\RI} j_{\ell}(k\RI)\right)  I_{3,\ell, m} \Bigg] .
  /// \f}
    /// The Clebsch-Gordan coefficients are calculated by calling the function c_ell_m(long, long) const. \n
  /// The time integrals \f$I_{\dots}\f$ are calculated when the function time_integration_theta_phi(long, long) is called.
  /// The function
  /// fill_bessel_array(double, double, double*) handles the computation of the Bessel functions \f$j_{\ell}(kr)\f$.
  void eval_spectrum_34_theta_phi() {
    if (labs(initial_m)>ell_grid_size) {
      cerr << "m is to large for N_{\ell}" << endl;
      exit(-1);
    };
    
    // in order to calculate the derivative of the spherical Bessel function j_l we use the j_{l+1} -> size of the array ell_grid_size+1
    double j_lkr[ell_grid_size+1];
    // big loop for summing over ells and such
    ofstream tsurff_polar_dat("tsurff-polar"+to_string(i_proc)+".dat");    
    tsurff_polar_dat.precision(17);
    for (long i_k=0; i_k<num_k_proc; i_k++) {
      const double k=k_values[i_k];
      fill_bessel_array(k, R_tsurff, j_lkr);
      for (long i_theta=0; i_theta<num_theta_surff; i_theta++) {
	const double theta=thetas_surff[i_theta];
	const double cos_theta=cos(thetas_surff[i_theta]);	
	cplxd b_k(0.0, 0.0);
	for (long i_ell=labs(initial_m); i_ell<ell_grid_size; i_ell++) {
	  // $c_{\ell, m}=\sqrt{\frac{(\ell+1)^2-m^2}{(2\ell+1)(2\ell+3)}}$
	  double c_ell_m=sqrt(pow2(i_ell+1)-pow2(initial_m))*sqrt(1./double((2*i_ell+1)*(2*i_ell+3)));
	  // \frac{d}{dr} f_{\ell}(kr)&=-kf_{\ell+1}(kr) + \frac{\ell}{r}f_{\ell}(kr))
	  const double j_lkr_deriv=-k*j_lkr[i_ell+1] + double(i_ell)/R_tsurff*j_lkr[i_ell];
	  const long ind_surff=index_ell_theta_k(i_ell, i_theta, i_k);

	  // spherical harmonics for \ell and \ell+1  
	  const cplxd Ylm=Ylm_prefactor(initial_m)*gsl_sf_legendre_sphPlm (i_ell, labs(initial_m), cos_theta);
	  const cplxd Ylp1m=Ylm_prefactor(initial_m)*gsl_sf_legendre_sphPlm (i_ell+1, labs(initial_m), cos_theta);
	  // contributions including the vector potential A
	  // \psi_{\ell} and Y_{\ell+1} are coupled here
	  const cplxd b_kl_Ylm_1=cplxd(0., 2.)*c_ell_m*pow_neg_i(i_ell+1)*j_lkr[i_ell+1]*Ylp1m*psi_A_surff[ind_surff];
	  // \psi_{\ell+1} and Y_{\ell} are coupled here (for \ell=N_{\ell}-1 this is set to zero)
	  const cplxd b_kl_psi_1=((i_ell+1)<ell_grid_size)?(cplxd(0., 2.)*c_ell_m*pow_neg_i(i_ell)*j_lkr[i_ell]*Ylm*psi_A_surff[ind_surff+1]):0.0;
	  const cplxd b_kl=b_kl_Ylm_1
	    +b_kl_psi_1
	    +pow_neg_i(i_ell)*j_lkr[i_ell]*Ylm*(-1.0/R_tsurff*psi_surff[ind_surff]+psi_deriv_surff[ind_surff])
	    -pow_neg_i(i_ell)*j_lkr_deriv*Ylm*psi_surff[ind_surff];
	  b_k+=b_kl;
	};
	b_k*=cplxd(0.0, -0.5)*4.*M_PI/sqrt(pow3(2.*M_PI))*R_tsurff*delta_t;
	// tsurff_polar_dat << k << " " << theta << " " << real(b_k) << " " << imag(b_k) << endl;
	tsurff_polar_dat << k*k*0.5 << " " << k << " " << theta << " " << norm(b_k)*k << endl;
      };
      // line break for every new theta
      tsurff_polar_dat << endl;
    };
  };

  /// calculate angle and energy resolved spectrum for polarization in \f$xy\f$-plane and expansion method 1 \n
  /// side effects: writes to file tsurff-polari_proc.dat \n
  /// reads time integrals from arrays psi_A_cc_surff, psi_A_surff psi_surff and psi_deriv_surff \n
  /// dependencies: Ylm_prefactor and gsl_sf_legendre_sphPlm for calculating the spherical harmonics \f$Y_{\ell m}\f$,
  ///               fill_bessel_array for calculating spherical bessel functions \n
  /// The probability amplitude is
  /// \f{equation*}{
  ///   \label{eq:pol-tsurff-ell-m-sum}
  /// b_1(\mathbf{k}) = \sum\limits^{N_{\ell}}_{\ell=0} \sum\limits^{\ell}_{m=-\ell} b_{1,\ell,m}(\mathbf{k}) Y_{\ell, m}(\Omega_k)
  /// \f}
  /// where
  /// \f{multline*}{
  ///   b_{1,\ell,m}(\mathbf{k})
  /// = \frac{1}{(2\pi)^{1/2}} \RI
  /// \Bigg[ \sqrt{2} {(-i)}^\ell j_\ell(k r) \\
  /// \times \Bigg\lbrace 
  /// - \left( I_{0,\ell-1, m-1} d_{--,\ell,m} - I_{0,\ell+1, m-1}  d_{+-,\ell,m} \right) 
  /// + \left( I_{1,\ell-1, m+1} d_{-+,\ell,m} - I_{1,\ell+1, m+1} d_{++,\ell,m} \right) \Bigg\rbrace \\
  /// +  (-i)^{\ell+1} j_{\ell}(k \RI)  
  /// \left( -\RI^{-1} I_{3,\ell,m} +  I_{4,\ell,m} \right)
  /// + (-i)^{\ell+1} k\left(j_{\ell+1}(k\RI) - \frac{\ell}{k\RI}j_{\ell}(k\RI)\right)  I_{3,\ell,m} \Bigg] .
  /// \f}
    /// The Clebsch-Gordan coefficients \f$d_{--},d_{-+},d_{+-},d_{++}\f$ are calculated by calling the functions
  /// clebsch_mm1_ellm1(long, long), clebsch_mp1_ellm1(long, long), clebsch_mm1_ellp1(long, long) and clebsch_mp1_ellp1(long, long). \n
  /// The time integrals \f$I_{\dots}\f$ are calculated when the function time_integration_theta_phi_ell_m(long, long) is called.
  /// The function
  /// fill_bessel_array(double, double, double*) handles the computation of the Bessel functions \f$j_{\ell}(kr)\f$.
  void eval_spectrum_44_theta_phi() {
    // in order to calculate the derivative of the spherical Bessel function j_l we use the j_{l+1} -> size of the array ell_grid_size+1
    double j_lkr[ell_grid_size+1];
    ofstream tsurff_polar_dat("tsurff-polar"+to_string(i_proc)+".dat");
    tsurff_polar_dat.precision(17);
    for (long i_k=0; i_k<num_k_proc; i_k++) {
      const double k=k_values[i_k];
      for (long i_theta=0; i_theta<num_theta_surff; i_theta++) {
	const double theta=thetas_surff[i_theta];
	const double cos_theta=cos(thetas_surff[i_theta]);
	for (long i_phi=0; i_phi<num_phi_surff; i_phi++) {
	  const double phi=phis_surff[i_phi];
	  // calculate spherical Bessel functions for fixed argument and all necessary \ell
	  fill_bessel_array(k, R_tsurff, j_lkr);
	  cplxd b_k(0.0, 0.0);
	  const double sqrt_2=sqrt(2.0);
	  for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	    // \frac{d}{dr} f_{\ell}(kr)&=-kf_{\ell+1}(kr) + \frac{\ell}{r}f_{\ell}(kr))
	    const double j_lkr_deriv=-k*j_lkr[i_ell+1] + double(i_ell)/R_tsurff*j_lkr[i_ell];
	    for (long i_m=-long(i_ell); i_m<long(i_ell+1); i_m++) {
	      const long ind_surff = index_ell_m_theta_phi_k(i_ell, i_m, i_theta, i_phi, i_k);

	      // Y_lm_prefactor is 1 if i_m \geq 0 or (-1)^m for i_m<0, Y_l,-m=(-1)^m Y_l,m^* 
	      const cplxd Ylm=Ylm_prefactor(i_m)*gsl_sf_legendre_sphPlm (i_ell, labs(i_m), cos_theta)*exp(cplxd(0.0, double(i_m)*phi));

	      const long N_ell = ell_grid_size;
	      if (i_ell>0 && i_m>=-long(i_ell)+2) {
		const long ind_surff_mm1_ellm1=index_ell_m_theta_phi_k(i_ell-1,i_m-1,i_theta,i_phi,i_k);
		b_k-=sqrt_2*pow_neg_i(i_ell)*j_lkr[i_ell]*Ylm*psi_A_cc_surff[ind_surff_mm1_ellm1]*clebsch_mm1_ellm1(i_ell, i_m);
	      };
	      if (i_ell>0 && i_m<=long(i_ell)-2) {
		const long ind_surff_mp1_ellm1=index_ell_m_theta_phi_k(i_ell-1,i_m+1,i_theta,i_phi,i_k);
		b_k+=sqrt_2*pow_neg_i(i_ell)*j_lkr[i_ell]*Ylm*psi_A_surff[ind_surff_mp1_ellm1]*clebsch_mp1_ellm1(i_ell, i_m);
	      };
	      if (i_ell<(N_ell-1)) {
		const long ind_surff_mm1_ellp1=index_ell_m_theta_phi_k(i_ell+1,i_m-1,i_theta,i_phi,i_k);
		b_k+=sqrt_2*pow_neg_i(i_ell)*j_lkr[i_ell]*Ylm*psi_A_cc_surff[ind_surff_mm1_ellp1]*clebsch_mm1_ellp1(i_ell, i_m);
	      };
	      if (i_ell<(N_ell-1)) {
		const long ind_surff_mp1_ellp1=index_ell_m_theta_phi_k(i_ell+1,i_m+1,i_theta,i_phi,i_k);
		b_k-=sqrt_2*pow_neg_i(i_ell)*j_lkr[i_ell]*Ylm*psi_A_surff[ind_surff_mp1_ellp1]*clebsch_mp1_ellp1(i_ell, i_m);
	      };
	      b_k+=pow_neg_i(i_ell+1)*j_lkr[i_ell]*Ylm*(-1.0/R_tsurff*psi_surff[ind_surff]+psi_deriv_surff[ind_surff])
		-pow_neg_i(i_ell+1)*j_lkr_deriv*Ylm*psi_surff[ind_surff];
	    };
	  };
	  b_k*=1.0/sqrt(2.*M_PI)*R_tsurff*delta_t;
	  // tsurff_polar_dat  << k << " " << theta << " " << phi << " " << real(b_k) << " " << imag(b_k) << endl;
	  tsurff_polar_dat << k*k*0.5 << " " << k << " " << theta << " " << phi << " " << norm(b_k)*k << endl;
	};
	// line break for every new phi
	tsurff_polar_dat << endl;
      };
    };
  };    
  
  inline long index_ell1_ell2(long i_ell1, long i_ell2) const {
    return  i_ell1 + i_ell2*ell_grid_size;
  };
  
  /// get the Wigner 3j coefficients and calculate necessary prefactors without any tricks (for linear polarization and fixed m)
  /// \f{multline*}
  /// \text{prefactor} \leftarrow w_{\ell,m,\ell_1,\ell_2} = \sqrt{(2\ell+1)(2\ell_1+1)} (2\ell_2+1) \\
  /// \times \begin{pmatrix}
  ///   \ell & \ell_1 & \ell_2\\
  ///   0 & 0 & 0
  /// \end{pmatrix}
  /// \begin{pmatrix}
  ///   \ell & \ell_1 & \ell_2\\
  ///   -m & m & 0
  /// \end{pmatrix} \quad \text{for all } \ell_1,\ell_2
  /// \f}
  void fill_wigner_prefactors(long ell, long m, double_ptr& prefactor) {
    for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
      for (long i_ell1=0; i_ell1<ell_grid_size; i_ell1++) {
	const double wig0 = gsl_sf_coupling_3j (2*ell, 2*i_ell1, 2*i_ell2, 0, 0, 0);
	const double wig1 = gsl_sf_coupling_3j (2*ell, 2*i_ell1, 2*i_ell2, -2*m, 2*m, 0);
	prefactor[index_ell1_ell2(i_ell1, i_ell2)]=sqrt(double((2*ell+1)*(2*i_ell1+1)))*double(2*i_ell2+1)*wig0*wig1;
      };
    };
  };
  
  /// A Clebsch-Gordan coefficient
  /// \f$c_{\ell, m}=\sqrt{\frac{(\ell+1)^2-m^2}{(2\ell+1)(2\ell+3)}}\f$
  inline double c_ell_m(long ell, long m) const {
#ifdef ADDITIONAL_TESTS
    if ((pow2(ell+1)-pow2(m))<0) {
      cerr << " c_ell_m " << ell << " " << m << endl;
      exit(-1);
    };
#endif
    return sqrt(pow2(ell+1)-pow2(m))*sqrt(1./double((2*ell+1)*(2*ell+3)));
  };

  inline long index_ell_k(long i_ell, long i_k) const {
    return  i_ell + i_k*ell_grid_size;
  };

  /// calculate the partial spectra for linear polarization and fixed \f$m\f$ \n
  /// side effects: writes result to array partial_amplitude \n
  /// reads time integrals from arrays psi_A_surff psi_surff and psi_deriv_surff
  /// \f{multline*}{
  /// \text{partial\_amplitude} \leftarrow b_{2,\ell, m}(k)= \frac{\RI}{(2\pi)^{1/2}} \sum_{\ell_1,\ell_2}
  /// {(-1)}^{m}
  /// {-i}^{\ell_1-\ell_2+1} w_{\ell,m,\ell_1,\ell_2} \\
  /// \times \Bigg[ 2 i j_{\ell_1}(k r)
  /// \left( c_{\ell_1-1,m} I_{2,\ell_1-1,m,\ell_2} + c_{\ell_1,m} I_{2,\ell_1+1,m,\ell_2} \right) \\
  /// + j_{\ell_1}(k r)
  /// \left(  I_{4,\ell_1,m,\ell_2} - \frac{1}{r} I_{3,\ell_1,m,\ell_2} \right)
  /// - I_{3,\ell_1,m,\ell_2} \partial_r j_{\ell_1}(k r) \Bigg]_{r=\RI}
  /// \f}
  /// The Clebsch-Gordan coefficients are calculated by calling the function c_ell_m(long, long) const. \n
  /// The time integrals \f$I_{\dots}\f$ are calculated when the function time_integration_ell1_ell2(long, long) is called.
  /// By calling the function fill_wigner_prefactors(long i_ell, long i_m, double_ptr &prefactor) the coefficients \f$w_{\ell,m,\ell_1,\ell_2}\f$ are computed and the function
  /// fill_bessel_array(double, double, double*) handles the computation of the Bessel functions \f$j_{\ell}(kr)\f$.
  void eval_partial_spectra_34_ell1_ell2() {
    // in order to calculate the derivative of the spherical Bessel function j_l we use the j_{l+1} -> size of the array ell_grid_size+1
    double j_lkr[ell_grid_size+1];
    double_ptr wigner_pre(new double[ell_grid_size*ell_grid_size]);
    // for every \ell and k calculate \sum_{\ell_1,\ell_2} w_{\ell,\ell_1,\ell_2,m} (I_0, ...)
    for (long i_ell=labs(initial_m); i_ell<ell_grid_size; i_ell++) {
      fill_wigner_prefactors(i_ell, initial_m, wigner_pre);
      for (long i_k=0; i_k<num_k_proc; i_k++) {
	const double k=k_values[i_k];
	fill_bessel_array(k, R_tsurff, j_lkr);
	cplxd b_kl_sum(0.0);
	for (long i_ell1=labs(initial_m); i_ell1<ell_grid_size; i_ell1++) {
	  // $c_{\ell, m}=\sqrt{\frac{(\ell+1)^2-m^2}{(2\ell+1)(2\ell+3)}}$
	  double c_ell1_m=c_ell_m(i_ell1, initial_m);
	  double c_ell1m1_m=c_ell_m(i_ell1-1, initial_m);
	  // \frac{d}{dr} f_{\ell}(kr)&=-kf_{\ell+1}(kr) + \frac{\ell}{r}f_{\ell}(kr))
	  const double j_lkr_deriv=-k*j_lkr[i_ell1+1] + double(i_ell1)/R_tsurff*j_lkr[i_ell1];
	  for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
	    const long ind_surff=index_ell1_ell2_k(i_ell1, i_ell2, i_k);
	    const long ind_surff_next=index_ell1_ell2_k(i_ell1+1, i_ell2, i_k);
	    const long ind_surff_prev=index_ell1_ell2_k(i_ell1-1, i_ell2, i_k);
	    const double w_ell1_ell2=wigner_pre[index_ell1_ell2(i_ell1, i_ell2)];
	    if (w_ell1_ell2!=0) {
	      b_kl_sum+=w_ell1_ell2*pow_neg_1(initial_m)*pow_neg_i(i_ell1+1)*pow_i(i_ell2)*
		(cplxd(0, 2.0)*j_lkr[i_ell1]*
		 (((i_ell1>labs(initial_m))?c_ell1m1_m*psi_A_surff[ind_surff_prev]:0.0) // if \ell_1 < m this part has to be zero
		  + (((i_ell1+1)<ell_grid_size)?c_ell1_m*psi_A_surff[ind_surff_next]:0.0)) // if \ell_1 = N_{\ell}-1 there is no \ell+1 orbital
		 + j_lkr[i_ell1]*(psi_deriv_surff[ind_surff]-psi_surff[ind_surff]/R_tsurff)
		 - j_lkr_deriv*psi_surff[ind_surff]);
	    };
	  };
	};
	partial_amplitude[index_ell_k(i_ell, i_k)] = b_kl_sum*R_tsurff*delta_t/sqrt(2.*M_PI);
      };
    };
  };

  /// calculate angle and energy resolved spectrum for linear polarization and fixed $m$
  /// side effects: writes to file tsurff-polari_proc.dat
  /// reads partial spectra from array partial_amplitude
  /// dependencies: Ylm_prefactor and gsl_sf_legendre_sphPlm for calculating the spherical harmonics \f$Y_{\ell m}\f$
  /// \f{equation*}{
  ///   b_2(\veck) = \sum_{\ell} b_{2,\ell, m}(k) Y_{\ell,m}(\Omega_k) 
  /// \f}  
  /// the coefficoents \f$b_{2,\ell, m}(k)\f$ are calculated when eval_partial_spectra_34_ell1_ell2() is called.
  void eval_spectrum_34_ell1_ell2() {
    ofstream tsurff_polar_dat("tsurff-polar"+to_string(i_proc)+".dat");    
    tsurff_polar_dat.precision(17);
    for (long i_k=0; i_k<num_k_proc; i_k++) {
      const double k=k_values[i_k];
      for (long i_theta=0; i_theta<num_theta_surff; i_theta++) {
	const double theta=thetas_surff[i_theta];
	const double cos_theta=cos(thetas_surff[i_theta]);
	cplxd b_k(0.0, 0.0);
	for (long i_ell=labs(initial_m); i_ell<ell_grid_size; i_ell++) {
	  const cplxd Ylm=Ylm_prefactor(initial_m)*gsl_sf_legendre_sphPlm (i_ell, labs(initial_m), cos_theta);
	  b_k+=partial_amplitude[index_ell_k(i_ell, i_k)]*Ylm;
	};
	// tsurff_polar_dat << k << " " << theta << " " << real(b_k) << " " << imag(b_k) << endl;
	tsurff_polar_dat << k*k*0.5 << " " << k << " " << theta << " " << norm(b_k)*k << endl;
      };
      // line break for every new theta
      tsurff_polar_dat << endl;
    };
  };

  inline long index_ell_m_k(long i_ell, long i_m, long i_k) const {
    const long ell_grid_size_2=ell_grid_size*ell_grid_size;
    return ((i_ell+1)*(i_ell+1)-(i_ell+1)+i_m) + i_k*ell_grid_size_2;
  };

  inline long index_ell1_m1_ell2_m2(long i_ell1, long i_m1, long i_ell2, long i_m2) const {
    const long ell_grid_size_2=ell_grid_size*ell_grid_size;
    return  ((i_ell1+1)*(i_ell1+1)-(i_ell1+1)+i_m1) + ((i_ell2+1)*(i_ell2+1)-(i_ell2+1)+i_m2)*ell_grid_size_2;
  };
  
  // calculate the prefactors without any tricks (inefficient)
  void fill_wigner_prefactors_ineff_44(long i_ell, long i_m, double_ptr& prefactor) {
    for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
      for (long i_m2=-(i_ell2); i_m2<(i_ell2+1); i_m2++) {
	for (long i_ell1=0; i_ell1<ell_grid_size; i_ell1++) {
	  for (long i_m1=-(i_ell1); i_m1<(i_ell1+1); i_m1++) {
	    const double wig0 = gsl_sf_coupling_3j (2*i_ell, 2*i_ell1, 2*i_ell2, 0, 0, 0);
	    const double wig1 = gsl_sf_coupling_3j (2*i_ell, 2*i_ell1, 2*i_ell2, -2*i_m, 2*i_m1, -2*i_m2);
	    prefactor[index_ell1_m1_ell2_m2(i_ell1, i_m1, i_ell2, i_m2)]=sqrt(double((2*i_ell+1)*(2*i_ell1+1)*(2*i_ell2+1)))*wig0*wig1;
	  };
	};
      };
    };
  };

  inline long index_ell1_m1_ell2(long i_ell1, long i_m1, long i_ell2) const {
    const long ell_grid_size_2=ell_grid_size*ell_grid_size;
    return  ((i_ell1+1)*(i_ell1+1)-(i_ell1+1)+i_m1) + i_ell2*ell_grid_size_2;
  };
  
  /// get the Wigner 3j coefficients and  calculate the necessary prefactors using -m+m_1-m_2=0 (for polarization in \f$xy\f$-plane)
  ///
  /// \f{multline*}
  /// \text{prefactor} \leftarrow w_{\ell,m,\ell_1,m_1,\ell_2,m_1-m} = \sqrt{(2\ell+1)(2\ell_1+1)(2\ell_2+1)} \\
  /// \times \begin{pmatrix}
  ///   \ell & \ell_1 & \ell_2\\
  ///   0 & 0 & 0
  /// \end{pmatrix}
  /// \begin{pmatrix}
  ///   \ell & \ell_1 & \ell_2\\
  ///   -m & m_1 & -(m_1-m)
  /// \end{pmatrix} \quad \text{for all } \ell_1,m_1,\ell_2
  /// \f}
  void fill_wigner_prefactors_44(long i_ell, long i_m, double_ptr& prefactor) {
    for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
      for (long i_ell1=0; i_ell1<ell_grid_size; i_ell1++) {
	for (long i_m1=-(i_ell1); i_m1<(i_ell1+1); i_m1++) {
	  const double wig0 = gsl_sf_coupling_3j (2*i_ell, 2*i_ell1, 2*i_ell2, 0, 0, 0);
	  const double wig1 = gsl_sf_coupling_3j (2*i_ell, 2*i_ell1, 2*i_ell2, -2*i_m, 2*i_m1, -2*(i_m1-i_m));
	  prefactor[index_ell1_m1_ell2(i_ell1, i_m1, i_ell2)]=sqrt(double((2*i_ell+1)*(2*i_ell1+1)*(2*i_ell2+1)))*wig0*wig1;
	};
      };
    };
  };

#define TESTS
  /// A Clebsch-Gordan coefficient
  ///
  /// \f$d_{--,\ell,m}=\sqrt{\frac{(\ell + m - 1)(\ell + m)}{2(2\ell - 1)(2\ell + 1)}}\f$
  inline double clebsch_mm1_ellm1(long ell, long m) {
    const double numerator=(ell + m - 1.0)*(ell + m);
    const double denominator=(2.0*(2.0*ell - 1.0)*(2.0*ell + 1.0));
#ifdef TESTS
    if (numerator/denominator < 0.0) {
      cerr << "negative sqrt in clebsch_mm1_ellm1 for m=" << m << " ell=" << ell << endl;      
    };
#endif
    return sqrt(numerator/denominator);
  };

  /// A Clebsch-Gordan coefficient
  ///
  /// \f$d_{+-,\ell,m}=\sqrt{\frac{(\ell - m + 1)(\ell - m + 2)(\ell + 1)}{(2\ell + 1)(2\ell + 2)(2\ell + 3)}}\f$
  inline double clebsch_mm1_ellp1(long ell, long m) {
    const double numerator=(ell - m + 1.0)*(ell - m + 2.0)*(ell + 1.0);
    const double denominator=(1.0*(2.0*ell + 1.0)*(2.0*ell + 2.0)*(2.0*ell + 3.0));
#ifdef TESTS
    if (numerator/denominator < 0.0) {
      cerr << "negative sqrt in clebsch_mm1_ellp1 for m=" << m << " ell=" << ell << endl;      
    };
#endif
    return sqrt(numerator/denominator);
  };

  /// A Clebsch-Gordan coefficient
  ///
  /// \f$d_{-+,\ell,m}=\sqrt{\frac{(\ell - m - 1)(\ell - m)}{2(2\ell - 1)(2\ell + 1)}}\f$
  inline double clebsch_mp1_ellm1(long ell, long m) {
    const double numerator=(ell - m - 1.0)*(ell - m);
    const double denominator=(2.0*(2.0*ell - 1.0)*(2.0*ell + 1.0));
#ifdef TESTS
    if (numerator/denominator < 0.0) {
      cerr << "negative sqrt in clebsch_mp1_ellm1 for m=" << m << " ell=" << ell << endl;      
    };
#endif
    return sqrt(numerator/denominator);
  };

  /// A Clebsch-Gordan coefficient
  ///
  /// \f$d_{++,\ell,m}=\sqrt{\frac{(\ell + m + 1)(\ell + m + 2)(\ell + 1)}{(2\ell + 1)(2\ell + 2)(2\ell + 3)}}\f$
  inline double clebsch_mp1_ellp1(long ell, long m) {
    const double numerator=(ell + m + 1.0)*(ell + m + 2.0)*(ell + 1.0);
    const double denominator=(1.0*(2.0*ell + 1.0)*(2.0*ell + 2.0)*(2.0*ell + 3.0));
#ifdef TESTS
    if (numerator/denominator < 0.0) {
      cerr << "negative sqrt in clebsch_mp1_ellp1 for m=" << m << " ell=" << ell << endl;
    };
#endif
    return sqrt(numerator/denominator);
  };

  /// calculate partial spectra for polarization in \f$xy\f$-plane and expansion method 2 \n
  /// side effects: write in array partial_amplitude \n
  /// reads time integrals from arrays psi_A_cc_surff, psi_A_surff psi_surff and psi_deriv_surff
  /// \f{multline*}{
  /// \text{partial\_amplitude} \leftarrow b_{2,\ell,m}(k) = \sqrt{2} \RI \sum_{\ell_1,\ell_2,m_1} w_{\ell,m,\ell_1,m_1,\ell_2,m_1-m} \\
  /// \times
  /// \Bigg\lbrace \sqrt{2} j_{\ell_1}(k r) 
  /// \Big[
  /// \left( I_{1,\ell_1-1,m_1+1,\ell_2,m_1-m} d_{-+,\ell_1,m_1} - I_{1,\ell_1+1,m_1+1,\ell_2,m_1-m} d_{++,\ell_1,m_1} \right) \\
  /// - \left(I_{0,\ell_1-1,m_1-1,\ell_2,m_1-m} d_{--,\ell_1,m_1}  - I_{0,\ell_1+1,m_1-1,\ell_2,m_1-m} d_{+-,\ell_1,m_1} \right)
  /// \Big] \\
  /// -i \left[ 
  /// j_{\ell_1}(k r) 
  /// \left(  I_{4,\ell_1,m_1,\ell_2,m_1-m} - \frac{1}{\RI} I_{3,\ell_1,m_1,\ell_2,m_1-m} \right) 
  /// -  I_{3,\ell_1,m_1,\ell_2,m_1-m} \partial_r j_{\ell_1}(k r) \right] \Bigg\rbrace_{r=\RI}
  /// \f}
  /// The Clebsch-Gordan coefficients \f$d_{--},d_{-+},d_{+-},d_{++}\f$ are calculated by calling the functions
  /// clebsch_mm1_ellm1(long, long), clebsch_mp1_ellm1(long, long), clebsch_mm1_ellp1(long, long) and clebsch_mp1_ellp1(long, long). \n
  /// The time integrals \f$I_{\dots}\f$ are calculated when the function time_integration_ell1_m1_ell2_m2(long, long) is called.
  /// By calling the function fill_wigner_prefactors_44(long i_ell, long i_m, double_ptr &prefactor) the coefficients \f$w_{\ell,m,\ell_1,m_1,\ell_2}\f$ are computed and the function
  /// fill_bessel_array(double, double, double*) handles the computation of the Bessel functions \f$j_{\ell}(kr)\f$.
  void eval_partial_spectra_44_ell1_ell2() {
    // in order to calculate the derivative of the spherical Bessel function j_l we use the j_{l+1} -> size of the array ell_grid_size+1
    double j_lkr[ell_grid_size+1];
    // for storage of coefficients for fixed \ell and m (one value for every combination of \ell_1, m_1 and \ell_2)
    double_ptr wigner_pre(new double[ell_grid_size*ell_grid_size*ell_grid_size]);
    // 
    for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
      for (long i_m=-(i_ell); i_m<(i_ell+1); i_m++) {
	fill_wigner_prefactors_44(i_ell, i_m, wigner_pre);
	for (long i_k=0; i_k<num_k_proc; i_k++) {
	  const double k=k_values[i_k];
	  fill_bessel_array(k, R_tsurff, j_lkr);
	  cplxd b_klm_sum(0.0);
	  for (long i_ell1=0; i_ell1<ell_grid_size; i_ell1++) {
	    for (long i_m1=-(i_ell1); i_m1<(i_ell1+1); i_m1++) {
	      // use selection rule for Wigner symbols
	      const long i_m2=i_m1-i_m;
	      // \frac{d}{dr} f_{\ell}(kr)&=-kf_{\ell+1}(kr) + \frac{\ell}{r}f_{\ell}(kr))
	      const double j_lkr_deriv=-k*j_lkr[i_ell1+1] + double(i_ell1)/R_tsurff*j_lkr[i_ell1];
	      for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
		const long ind_surff=index_ell1_m1_ell2_m2_k(i_ell1, i_m1, i_ell2, i_m2, i_k);
		const double w_ellmell1m1ell2=wigner_pre[index_ell1_m1_ell2(i_ell1, i_m1, i_ell2)];
		if (w_ellmell1m1ell2!=0) {
		  // cout << i_ell << " " << i_m << " " << i_ell1 << " " << i_m1 << " " << i_ell2 << " " << i_m2 << " " << w_ellmell1m1ell2 << endl;
		  const long N_ell = ell_grid_size;
		  if (i_ell1>0 && i_m1>=-long(i_ell1)+2) {
		    const long ind_surff_mm1_ellm1=index_ell1_m1_ell2_m2_k(i_m1-1, long(i_ell1)-1, i_ell2, i_m2, i_k);
		    b_klm_sum-=w_ellmell1m1ell2*pow_neg_1(i_m1+i_m2)*pow_i(i_ell2-i_ell1)*j_lkr[i_ell1]*psi_A_cc_surff[ind_surff_mm1_ellm1]*clebsch_mm1_ellm1(i_ell1, i_m1);
		  };
		  if (i_ell1>0 && i_m1<=long(i_ell1)-2) {
		    const long ind_surff_mp1_ellm1=index_ell1_m1_ell2_m2_k(i_m1+1, long(i_ell1)-1, i_ell2, i_m2, i_k);
		    b_klm_sum+=w_ellmell1m1ell2*pow_neg_1(i_m1+i_m2)*pow_i(i_ell2-i_ell1)*j_lkr[i_ell1]*psi_A_surff[ind_surff_mp1_ellm1]*clebsch_mp1_ellm1(i_ell1, i_m1);
		  };
		  if (i_ell1<(N_ell-1)) {
		    const long ind_surff_mm1_ellp1=index_ell1_m1_ell2_m2_k(i_m1-1, i_ell1+1, i_ell2, i_m2, i_k);
		    b_klm_sum+=w_ellmell1m1ell2*pow_neg_1(i_m1+i_m2)*pow_i(i_ell2-i_ell1)*j_lkr[i_ell1]*psi_A_cc_surff[ind_surff_mm1_ellp1]*clebsch_mm1_ellp1(i_ell1, i_m1);
		  };
		  if (i_ell1<(N_ell-1)) {
		    const long ind_surff_mp1_ellp1=index_ell1_m1_ell2_m2_k(i_m1+1, i_ell1+1, i_ell2, i_m2, i_k);
		    b_klm_sum-=w_ellmell1m1ell2*pow_neg_1(i_m1+i_m2)*pow_i(i_ell2-i_ell1)*j_lkr[i_ell1]*psi_A_surff[ind_surff_mp1_ellp1]*clebsch_mp1_ellp1(i_ell1, i_m1);
		  };
		
		  b_klm_sum+=cplxd(0,-sqrt(0.5))*w_ellmell1m1ell2*pow_neg_1(i_m1+i_m2)*pow_i(i_ell2-i_ell1)*
		    (j_lkr[i_ell1]*(-1.0/R_tsurff*psi_surff[ind_surff] + psi_deriv_surff[ind_surff])-j_lkr_deriv*psi_surff[ind_surff]);
		};
	      };
	    };
	  };
	  partial_amplitude[index_ell_m_k(i_ell, i_m, i_k)] = b_klm_sum*2.0*R_tsurff*delta_t;
	};
      };
    };
  };

  /// calculate angle and energy resolved spectrum for polarization in \f$xy\f$-plane \n
  /// side effects: writes to file tsurff-polari_proc.dat \n
  /// reads partial spectra from array ::partial_amplitude \n
  /// dependencies: Ylm_prefactor and gsl_sf_legendre_sphPlm for calculating the spherical harmonics \f$Y_{\ell m}\f$
  /// \f{equation*}{
  ///   b_2(\veck)=\sum_{\ell, m} b_{2,\ell, m}(k) Y_{\ell,m}(\Omega_k)
  /// \f}
  /// the coefficoents \f$b_{2,\ell, m}(k)\f$ are calculated when eval_partial_spectra_44_ell1_ell2() is called.
  void eval_spectrum_44_ell1_ell2() {
    // big ass loop for summing over ells and such
    ofstream tsurff_polar_dat("tsurff-polar"+to_string(i_proc)+".dat");    
    tsurff_polar_dat.precision(17);
    for (long i_k=0; i_k<num_k_proc; i_k++) {
      const double k=k_values[i_k];
      for (long i_theta=0; i_theta<num_theta_surff; i_theta++) {
	const double theta=thetas_surff[i_theta];
	const double cos_theta=cos(thetas_surff[i_theta]);
	for (long i_phi=0; i_phi<num_phi_surff; i_phi++) {
	  const double phi=phis_surff[i_phi];
	  cplxd b_k(0.0, 0.0);
	  for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	    for (long i_m=-(i_ell); i_m<(i_ell+1); i_m++) {
	      const cplxd Ylm=Ylm_prefactor(i_m)*gsl_sf_legendre_sphPlm (i_ell, labs(i_m), cos_theta)*exp(cplxd(0.0, double(i_m)*phi));
	      b_k+=partial_amplitude[index_ell_m_k(i_ell, i_m, i_k)]*Ylm;
	    };
	  };
	  // tsurff_polar_dat << k*k*0.5 << " " << k << " " << theta << " " << phi << " " << real(b_k) << " " << imag(b_k) << endl;
	  tsurff_polar_dat << k*k*0.5 << " " << k << " " << theta << " " << phi << " " << norm(b_k)*k << endl;
	};
	// line break for every new phi
	tsurff_polar_dat << endl;
      };
    };
  };

  // functions for debugging purposes
public:
  int test_index_ell1_m1_ell2() {
    long counter;
    for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
      for (long i_ell1=0; i_ell1<ell_grid_size; i_ell1++) {
	for (long i_m1=-(i_ell1); i_m1<(i_ell1+1); i_m1++) {
	  if (counter!=index_ell1_m1_ell2(i_ell1, i_m1, i_ell2))
	    return counter;
	  counter++;
	};
      };
    };
    return 0;
  };


  /// for debugging purposes: write the Legendre polynomials to a file
  void print_Plms() {
    if (i_proc==0) {
      ofstream debug_Plm_dat("debug-Plm.dat");
      long num_Plm_one_theta=gsl_sf_legendre_array_size(ell_grid_size, initial_m);
      debug_Plm_dat << num_Plm_one_theta << " num_Plm_one_theta" << endl;
      double Plm_theta[num_Plm_one_theta*num_theta_surff];
      for (long i_theta=0; i_theta<num_theta_surff; i_theta++) {
	const double cos_theta=cos(thetas_surff[i_theta]);
	const double theta=(thetas_surff[i_theta]);
	// gsl_sf_legendre_sphPlm_array(ell_grid_size, initial_m, cos_theta, &(Plm_theta[i_theta*num_Plm_one_theta]));
	for (long i_ell=labs(initial_m); i_ell<ell_grid_size-1; i_ell++) {
	  // const long ind_Plm=i_ell-labs(initial_m)+num_Plm_one_theta*i_theta;
	  const double gsl_Plm = gsl_sf_legendre_sphPlm(i_ell, labs(initial_m), cos_theta);
	  // const cplxd dieter_Plm = ylm2(i_ell, initial_m, theta, 0.0);
	  debug_Plm_dat << i_theta << " " << i_ell << " " << gsl_Plm << endl;
	};
	debug_Plm_dat << endl;
      };
    };
  };

  /// For debugging purposes: print values of the Bessel function
  void print_bessel() {
    if (i_proc==0) {
      ofstream debug_bessel_dat("debug-bessel.dat");
      double j_l2kalpha[ell_grid_size];
      double j_l2kalpha_neg[ell_grid_size];
      for (long i_k=0; i_k<num_k_surff; i_k++) {
	const double k=k_values[i_k];
	fill_bessel_array(k, 1.0, j_l2kalpha);
	fill_bessel_array(k, -1.0, j_l2kalpha_neg);
	for (long i_ell_2=0; i_ell_2<ell_grid_size; i_ell_2++) {
	  debug_bessel_dat << i_ell_2 << " " << k << " " << j_l2kalpha[i_ell_2] << " " << j_l2kalpha_neg[i_ell_2] << endl;
	};
	debug_bessel_dat << endl;
      };
    };
  };
  
  /// For debugging purposes: print values of the time integrals psi\_surff, psi\_A\_surff, ...
  void print_int_dt_psi(long i_theta=0) {
    if (expansion_scheme==1) {
      ofstream tsurff_input_dat("tsurff-int_dt_psi-"+to_string(i_proc)+".dat");
      const double cos_theta=cos(thetas_surff[i_theta]);
      for (long i_k=0; i_k<num_k_proc; i_k++) {
	const double k=k_values[i_k];
	for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	  const long ind_surff=index_ell_theta_k(i_ell, i_theta, i_k);
	  tsurff_input_dat << i_ell << " " << k << " " << real(psi_surff[ind_surff]) << " " << imag(psi_surff[ind_surff]) << " "
			   << real(psi_deriv_surff[ind_surff]) << " " << imag(psi_deriv_surff[ind_surff]) << " "
			   << real(psi_A_surff[ind_surff]) << " " << imag(psi_A_surff[ind_surff]) << endl;
	};
	tsurff_input_dat << endl;
      };
    }
    else if (expansion_scheme==2 && qprop_dim==34) {
      ofstream tsurff_input_dat("tsurff-int_dt_psi"+to_string(i_proc)+".dat");
      for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
	tsurff_input_dat << "# ell_2=" << i_ell2 << endl;
	for (long i_k=0; i_k<num_k_proc; i_k++) {
	  const double k=k_values[i_k];
	  for (long i_ell1=0; i_ell1<ell_grid_size; i_ell1++) {
	    const long ind_surff=index_ell1_ell2_k(i_ell1, i_ell2, i_k);
	    tsurff_input_dat << i_ell1 << " " << k << " " << real(psi_surff[ind_surff]) << " " << imag(psi_surff[ind_surff]) << " "
			     << real(psi_deriv_surff[ind_surff]) << " " << imag(psi_deriv_surff[ind_surff]) << " "
			     << real(psi_A_surff[ind_surff]) << " " << imag(psi_A_surff[ind_surff]) << endl;
	  };
	tsurff_input_dat << endl;
	};
      tsurff_input_dat << endl;
      };
    }
    else if (expansion_scheme==2 && qprop_dim==44) {
      ofstream tsurff_input_dat("tsurff-int_dt_psi"+to_string(i_proc)+".dat");
      for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
	for (long i_m2=-(i_ell2); i_m2<(i_ell2+1); i_m2++) {
	  tsurff_input_dat << "# ell_2=" << i_ell2 << " m_2=" << i_m2 << endl;
	  for (long i_k=0; i_k<num_k_proc; i_k++) {
	    const double k=k_values[i_k];
	    for (long i_ell1=0; i_ell1<ell_grid_size; i_ell1++) {
	      tsurff_input_dat << i_ell1 << " " << k << " ";
	      for (long i_m1=-(i_ell2); i_m1<(i_ell2+1); i_m1++) {
		const long ind_surff=index_ell1_m1_ell2_m2_k(i_ell1, i_m1, i_ell2, i_m2, i_k);
		tsurff_input_dat << i_m1 << " " << real(psi_surff[ind_surff]) << " " << imag(psi_surff[ind_surff]) << " "
				 << real(psi_deriv_surff[ind_surff]) << " " << imag(psi_deriv_surff[ind_surff]) << " "
				 << real(psi_A_surff[ind_surff]) << " " << imag(psi_A_surff[ind_surff]) << " ";
	      };
	      tsurff_input_dat << endl;
	    };
	    tsurff_input_dat << endl;
	  };
	  tsurff_input_dat << endl;
	};
      };
    };
  };

  void print_wigner() {
    if (i_proc==0) {
      ofstream debug_wigner_dat("debug-wigner.dat");
      if (qprop_dim==34) {
	double_ptr wigner_pre(new double[ell_grid_size*ell_grid_size]);
	for (long i_ell=0; i_ell<ell_grid_size; i_ell++) {
	  debug_wigner_dat << "# ell=" << i_ell << endl;
	  fill_wigner_prefactors(i_ell, initial_m, wigner_pre);
	  for (long i_ell2=0; i_ell2<ell_grid_size; i_ell2++) {
	    for (long i_ell1=0; i_ell1<ell_grid_size; i_ell1++) {
	      debug_wigner_dat << i_ell2 << " " <<  i_ell1 << " " << wigner_pre[index_ell1_ell2(i_ell1, i_ell2)] << endl;
	    };
	    debug_wigner_dat << endl;
	  };
	  debug_wigner_dat << endl;
	};
      };
    };
  };
};

#endif
