#ifndef _VECPOT_H_
#define _VECPOT_H_

#include <iostream>
#include <cmath>

#include "powers.hh"
#include "grid.h"

double always_zero2(double t, int me);
double always_zero5(double x, double y, double z, double t, int me);
double always_zero_imag(long x, long y, long z, double t, grid g); 


// The vector potential with sine squared envelope
class vecpot {
  double omega;
  double n_c;
  double E_0;
  double phi_cep;
  double duration;
  double ww;
  double start_time, end_time;
public:
  vecpot() { duration = -1; };
  vecpot(double om, double n_cyc, double E_max, double cep, double start_time_in = 0.0) : omega(om), n_c(n_cyc), E_0(E_max), phi_cep(cep), start_time(start_time_in) {
    duration=n_c*2*M_PI/omega;
    // angular frequency of the envelope
    ww=omega/(2.0*n_c);
    end_time = start_time + duration;
  };
  double operator()(double time, int me) const {
    // switch of the field for propagation after the pulse
    if ((time>0.0) && (time<duration)) {
      return E_0/omega*pow2(sin(ww*time))*sin(omega*time+phi_cep); // here's the shape of the laser pulse
    }
    else {
      return 0;
    };
  };
  double get_omega() { return omega; }
  double get_E0() { return E_0; }
  double get_num_cycles() { return n_c; }
  double get_phase() { return phi_cep; }
  double get_start_time() { return start_time; }
  double get_end_time() { return end_time; }
  double get_duration() {
    return duration;
  };
  double get_Up() {
    return (E_0 * E_0)/4.0/(omega*omega);
  };
  double integral(double Time) {
    // grind(ratsimp(integrate(sin(pi*time/pulse_dur)**2*sin(omega*time+phi), time, 0, Time)));
    // const double omega = omega;
    const double omega_2 = pow2(omega);
    const double omega_3 = pow3(omega);
    const double pulse_dur = n_c*2.0*M_PI/omega;
    const double pulse_dur_2 = pow2(pulse_dur);
    const double pi_2 = pow2(M_PI);
    const double pi = M_PI;
    const double sin_phi_2 = pow2(sin(phi_cep));
    const double cos_phi_2 = pow2(cos(phi_cep));
    const double phi = phi_cep;
    const double ampl = E_0/omega;
    if (Time<pulse_dur) {
      return ampl*((omega_2*sin(phi)*pulse_dur_2-2.0*omega*sin(phi)*pi*pulse_dur)*sin(((omega*pulse_dur+2.0*pi)*Time+2.0*phi*pulse_dur)/pulse_dur)
		   +(omega_2*cos(phi)*pulse_dur_2-2.0*omega*cos(phi)*pi*pulse_dur)*cos(((omega*pulse_dur+2.0*pi)*Time+2.0*phi*pulse_dur)/pulse_dur)
		   +(omega_2*sin(phi)*pulse_dur_2+2.0*omega*sin(phi)*pi*pulse_dur)*sin(((omega*pulse_dur-2.0*pi)*Time+2.0*phi*pulse_dur)/pulse_dur)
		   +(omega_2*cos(phi)*pulse_dur_2+2.0*omega*cos(phi)*pi*pulse_dur)*cos(((omega*pulse_dur-2.0*pi)*Time+2.0*phi*pulse_dur)/pulse_dur)
		   +(8.0*sin(phi)*pi_2-2.0*omega_2*sin(phi)*pulse_dur_2)*sin(omega*Time+2.0*phi)
		   +(8.0*cos(phi)*pi_2-2.0*omega_2*cos(phi)*pulse_dur_2)*cos(omega*Time+2.0*phi)
		   +(2.0*omega*sin(phi)*pi*pulse_dur-omega_2*sin(phi)*pulse_dur_2)*sin((omega*pulse_dur+2.0*pi)*Time/pulse_dur)
		   +(omega_2*cos(phi)*pulse_dur_2-2.0*omega*cos(phi)*pi*pulse_dur)*cos((omega*pulse_dur+2.0*pi)*Time/pulse_dur)
		   +(-omega_2*sin(phi)*pulse_dur_2-2.0*omega*sin(phi)*pi*pulse_dur)*sin((omega*pulse_dur-2.0*pi)*Time/pulse_dur)
		   +(omega_2*cos(phi)*pulse_dur_2+2.0*omega*cos(phi)*pi*pulse_dur)*cos((omega*pulse_dur-2.0*pi)*Time/pulse_dur)
		   +(2.0*omega_2*sin(phi)*pulse_dur_2-8.0*sin(phi)*pi_2)*sin(omega*Time)
		   +(8.0*cos(phi)*pi_2-2.0*omega_2*cos(phi)*pulse_dur_2)*cos(omega*Time)
		   +(-8.0*sin(phi)*sin(2.0*phi)-8.0*cos(phi)*cos(2.0*phi)-8.0*cos(phi))*pi_2)
	/((8.0*omega_3*sin_phi_2+8.0*omega_3*cos_phi_2)*pulse_dur_2+(-32.0*omega*sin_phi_2-32.0*omega*cos_phi_2)*pi_2);
    }
    else {
      return ampl*((omega_2*sin(phi)*pulse_dur_2-2.0*omega*sin(phi)*pi*pulse_dur)*sin(((omega*pulse_dur+2.0*pi)*pulse_dur+2.0*phi*pulse_dur)/pulse_dur)
		   +(omega_2*cos(phi)*pulse_dur_2-2.0*omega*cos(phi)*pi*pulse_dur)*cos(((omega*pulse_dur+2.0*pi)*pulse_dur+2.0*phi*pulse_dur)/pulse_dur)
		   +(omega_2*sin(phi)*pulse_dur_2+2.0*omega*sin(phi)*pi*pulse_dur)*sin(((omega*pulse_dur-2.0*pi)*pulse_dur+2.0*phi*pulse_dur)/pulse_dur)
		   +(omega_2*cos(phi)*pulse_dur_2+2.0*omega*cos(phi)*pi*pulse_dur)*cos(((omega*pulse_dur-2.0*pi)*pulse_dur+2.0*phi*pulse_dur)/pulse_dur)
		   +(8.0*sin(phi)*pi_2-2.0*omega_2*sin(phi)*pulse_dur_2)*sin(omega*pulse_dur+2.0*phi)
		   +(8.0*cos(phi)*pi_2-2.0*omega_2*cos(phi)*pulse_dur_2)*cos(omega*pulse_dur+2.0*phi)
		   +(2.0*omega*sin(phi)*pi*pulse_dur-omega_2*sin(phi)*pulse_dur_2)*sin((omega*pulse_dur+2.0*pi)*pulse_dur/pulse_dur)
		   +(omega_2*cos(phi)*pulse_dur_2-2.0*omega*cos(phi)*pi*pulse_dur)*cos((omega*pulse_dur+2.0*pi)*pulse_dur/pulse_dur)
		   +(-omega_2*sin(phi)*pulse_dur_2-2.0*omega*sin(phi)*pi*pulse_dur)*sin((omega*pulse_dur-2.0*pi)*pulse_dur/pulse_dur)
		   +(omega_2*cos(phi)*pulse_dur_2+2.0*omega*cos(phi)*pi*pulse_dur)*cos((omega*pulse_dur-2.0*pi)*pulse_dur/pulse_dur)
		   +(2.0*omega_2*sin(phi)*pulse_dur_2-8.0*sin(phi)*pi_2)*sin(omega*pulse_dur)
		   +(8.0*cos(phi)*pi_2-2.0*omega_2*cos(phi)*pulse_dur_2)*cos(omega*pulse_dur)
		   +(-8.0*sin(phi)*sin(2.0*phi)-8.0*cos(phi)*cos(2.0*phi)-8.0*cos(phi))*pi_2)
	/((8.0*omega_3*sin_phi_2+8.0*omega_3*cos_phi_2)*pulse_dur_2+(-32.0*omega*sin_phi_2-32.0*omega*cos_phi_2)*pi_2); // this is the brobdingnagian time integral of A(t), i.e., the free-electron excursion alpha(t)
    };
  };
};                                                                                   



//class scalarpot {
//  double nuclear_charge;
//  double R_co;
//public:
//  scalarpot(double charge, double co) : nuclear_charge(charge), R_co(co) {
//    //
//  };
//  double operator()(double x, double y, double z, double time, int me) const {
//    // scrinzi potential
//    // double result=(x<R_co)?nuclear_charge*(-1.0/x-pow2(x)/(2*pow3(R_co))+3.0/(2.0*R_co)):0.0;
//    // stupid simple Volker potential; first -1/r then linear then zero
//    const double m=1.0/pow2(R_co);
//    double result=(x<R_co)?-1.0/x:((x<2*R_co)?-1.0/R_co+m*(x-R_co):0.0);
//    return result;
//  };
//  double get_nuclear_charge() {return nuclear_charge; };
//};


// class scalarpot_alpha : public scalarpot {
class scalarpot {

  private:
  double nuclear_charge;
  double R_co;
  double alpha;

  public:
  scalarpot(double charge, double co) :
    nuclear_charge(charge), R_co(co), alpha(1.0) {}

  public:
  scalarpot(double charge, double co, double alpha_in) : 
    nuclear_charge(charge), R_co(co), alpha(alpha_in) {}

  double get_nuclear_charge() { return nuclear_charge; };

  double operator () (double x, double y, double z, double time, int me) const {
		double Z = nuclear_charge;
    double result = 0;
		double result0 = - (1 + exp(-alpha*R_co) * (Z - 1)) / R_co;
    if (x < R_co) {
			result = - (1 + exp(-alpha*x) * (Z - 1)) / x;
    } else if (x < 2*R_co) {
      result = - (x - R_co) / (R_co)*result0 + result0;
    } else {
      result = 0;
    }
    return result;
  }
};


class imagpot {
  double ampl_im;  // amplitude of imaginary absorbing potential  <--------------  100.0 imag pot on,  0.0 off
  long imag_potential_width;
public:
  imagpot(long width, double ampl=100.0) : ampl_im(ampl), imag_potential_width(width) {
    //
  };
  double operator()(long xindex, long yindex, long zindex, double time, grid g) {
    return (*this)(xindex, yindex, zindex, time, g.ngps_x());
  }
  double operator()(long xindex, long yindex, long zindex, double time, long ngps_x) {
    if (ampl_im>1.0) {
      const long imag_start=ngps_x-imag_potential_width;
      if (xindex<imag_start)
	return 0;
      else {
	const double r=double(xindex-imag_start)/double(imag_potential_width);
	return ampl_im*pow2(pow2(pow2(pow2(r))));
      };
    }
    else {
      return 0.0;
    };
  };
};


class superposed_vecpot {
  private:
  vecpot *vecpot_arr;
  long num_of_vecpot;

  public:
  superposed_vecpot() {}
  superposed_vecpot(vecpot *vecpot_arr_in, long num_of_vecpot_in) {

    // Check input arguments
    if (num_of_vecpot_in < 1) { 
      std::cerr << "[ERROR] Unexpected number of vecpot: "
        << num_of_vecpot << endl; exit(-1); }

    // Assigen arguments into member variables
    vecpot_arr = vecpot_arr_in;
    num_of_vecpot = num_of_vecpot_in;
  }

  double operator () (double time, int me) {
    double result = 0;
    long vecpot_index;
    for (vecpot_index = 0; vecpot_index < num_of_vecpot; vecpot_index++) {
      result += vecpot_arr[vecpot_index](time, me);
    }
    return result;
  }

  double integral(double Time) {
    double result = 0;
    long vecpot_index;
    for (vecpot_index = 0; vecpot_index < num_of_vecpot; vecpot_index++) {
      result += vecpot_arr[vecpot_index].integral(Time);
    }
    return result;
  }

  double get_start_time() {
    vecpot vp = vecpot_arr[0];
    long vecpot_index;
    double min_start_time = vp.get_start_time();
    for (vecpot_index = 1; vecpot_index < num_of_vecpot; vecpot_index++) {
      vp = vecpot_arr[vecpot_index];
      if (min_start_time > vp.get_start_time()) { min_start_time = vp.get_start_time(); }
    }
    return min_start_time;
  }

  double get_end_time() {
    vecpot vp = vecpot_arr[0];
    long vecpot_index;
    double max_end_time = vp.get_end_time();
    for (vecpot_index = 1; vecpot_index < num_of_vecpot; vecpot_index++) {
      vp = vecpot_arr[vecpot_index];
      if (max_end_time < vp.get_end_time()) { max_end_time = vp.get_end_time(); }
    }
    return max_end_time;
  }

  double get_duration() { return get_end_time() - get_start_time(); }
  long get_num_of_vecpot() { return num_of_vecpot; }

  vecpot *get_vecpot_arr() { return vecpot_arr; }

  // double operator += () ...
};

#endif  // _VECPOT_H_
