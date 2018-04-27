#include <memory>
#include <complex>
#include <string>
#include <iostream>

typedef std::complex<double> cplxd;
typedef std::unique_ptr<cplxd[]> cplxd_ptr;

#include <tsurffSpectrum.hh>
#include <potentials.hh>

#include <vecpot.hh>

using std::cout;
using std::string;

int main(int argc, char **argv) {

  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");

  superposed_vecpot vecpot_x, vecpot_y, vecpot_z;
  construct_vecpot(para_ini.getLong("qprop-dim"), para_prop, vecpot_x, vecpot_y, vecpot_z);

//  const double n_cyc=para_prop.getDouble("num-cycles");
//  vecpot vecpot_x(para_prop.getDouble("omega"), n_cyc, 0.0, 0.0);
//  vecpot vecpot_y(para_prop.getDouble("omega"), n_cyc, 0.0, 0.0);
//  vecpot vecpot_z(para_prop.getDouble("omega"), n_cyc, para_prop.getDouble("max-electric-field"), 0.0);

  
  // [NOTE] This should be modified to use other type, other than vecpot -> superposed_vecpot etc.
//  tsurffSpectrum<vecpot, vecpot, vecpot> tsurff(para_ini, para_prop, para_tsurff,
//      vecpot_x, vecpot_y, vecpot_z);
  print_vecpot(vecpot_x.get_vecpot_arr()[0], "in tsurff, vecpot_x");
  print_vecpot(vecpot_z.get_vecpot_arr()[0], "in tsurff, vecpot_z");
  
  cout << "vpz duration: " << vecpot_z.get_duration() << endl;  

  tsurffSpectrum<superposed_vecpot, superposed_vecpot, superposed_vecpot>
    tsurff(para_ini, para_prop, para_tsurff, vecpot_x, vecpot_y, vecpot_z);

  tsurff.time_integration();
  tsurff.polar_spectrum();
  tsurff.print_partial_amplitudes();

  return 0;
};
