// definitions of potentials

double always_zero2(double t, int me) {
  return 0;
};

double always_zero5(double x, double y, double z, double t, int me) {
  return 0;
};

double always_zero_imag(long x, long y, long z, double t, grid g) {
  return 0;
};

// The vector potential with sine squared envelope
class vecpot {
  double n_c;
  double phi_cep;
  double omega;
  double E_0;
  double duration;
  double ww;
public:
  vecpot(double om, double n_cyc, double E_max, double cep) : omega(om), n_c(n_cyc), E_0(E_max), phi_cep(cep) {
    duration=n_c*2*M_PI/omega;
    // angular frequency of the envelope
    ww=omega/(2.0*n_c);
  };
  double operator()(double time, int me) const {
    // switch of the field for propagation after the pulse
    if ((time>0.0) && (time<duration)) {
      return E_0/omega*pow2(sin(ww*time))*sin(omega*time+phi_cep);
    }
    else {
      return 0;
    };
  };
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
	/((8.0*omega_3*sin_phi_2+8.0*omega_3*cos_phi_2)*pulse_dur_2+(-32.0*omega*sin_phi_2-32.0*omega*cos_phi_2)*pi_2);
    };
  };
};                                                                                   

class scalarpot {
  double nuclear_charge;
  double R_co;
public:
  scalarpot(double charge, double co) : nuclear_charge(charge), R_co(co) {
    //
  };
  double operator()(double x, double y, double z, double time, int me) const {
    // scrinzi potential
    // double result=(x<R_co)?nuclear_charge*(-1.0/x-pow2(x)/(2*pow3(R_co))+3.0/(2.0*R_co)):0.0;
    // stupid simple Volker potential; first -1/r then linear then zero
    const double m=1.0/pow2(R_co);
    double result=(x<R_co)?-1.0/x:((x<2*R_co)?-1.0/R_co+m*(x-R_co):0.0);
    return result;
  };
  double get_nuclear_charge() {return nuclear_charge; };
};

class imagpot {
  long imag_potential_width;
  double ampl_im;  // amplitude of imaginary absorbing potential  <--------------  100.0 imag pot on,  0.0 off
public:
  imagpot(long width, double ampl=100.0) : ampl_im(ampl), imag_potential_width(width) {
    //
  };
  double operator()(long xindex, long yindex, long zindex, double time, grid g) {
    if (ampl_im>1.0) {
      const long imag_start=g.ngps_x()-imag_potential_width;
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
